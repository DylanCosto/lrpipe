//! Reader for PyTorch `torch.save` checkpoints holding a state dict, as the original ships them.
//!
//! The file is a zip archive with `<name>/data.pkl` (a pickle of an `OrderedDict` mapping
//! parameter names to `torch._utils._rebuild_tensor_v2(...)` calls) and one raw little-endian
//! file per storage in `<name>/data/<key>`. Only the pickle subset `torch.save` emits for such a
//! dict is understood; anything else is an error rather than a guess.

use std::collections::HashMap;
use std::path::Path;

use crate::Error;

/// Element data of one checkpoint tensor, in its storage's element type.
#[derive(Debug, Clone, PartialEq)]
pub enum TensorValues {
    /// `torch.float32`.
    F32(Vec<f32>),
    /// `torch.int64`.
    I64(Vec<i64>),
}

/// One tensor of a state dict, made contiguous.
#[derive(Debug, Clone, PartialEq)]
pub struct CheckpointTensor {
    /// Sizes of each dimension.
    pub shape: Vec<usize>,
    /// Elements in row-major order.
    pub values: TensorValues,
}

/// The named tensors of a state-dict checkpoint, in file order.
///
/// # Errors
///
/// Fails if the file is not a zip archive, uses pickle features outside what `torch.save` writes
/// for a state dict, or references missing storages.
pub fn read_state_dict(path: &Path) -> Result<Vec<(String, CheckpointTensor)>, Error> {
    let bad = |what: String| Error::Checkpoint { path: path.to_owned(), what };
    let data = std::fs::read(path).map_err(|e| bad(e.to_string()))?;
    let zip = StoredZip::parse(&data).map_err(&bad)?;
    let pkl_name = zip
        .entries
        .keys()
        .filter(|n| n.ends_with("/data.pkl") || *n == "data.pkl")
        .min()
        .ok_or_else(|| bad("no data.pkl".into()))?;
    let prefix = pkl_name.strip_suffix("data.pkl").expect("matched above");
    let read = |name: &str| zip.get(&data, name).map_err(|what| bad(format!("{name}: {what}")));
    if let Ok(order) = read(&format!("{prefix}byteorder"))
        && order != b"little"
    {
        return Err(bad(format!("unsupported byte order {:?}", String::from_utf8_lossy(order))));
    }
    let entries = Unpickler::new(read(pkl_name)?).run().map_err(&bad)?;
    entries
        .into_iter()
        .map(|(name, spec)| {
            let bytes = read(&format!("{prefix}data/{}", spec.storage.key))?;
            let tensor = spec.materialize(bytes).map_err(|what| bad(format!("{name}: {what}")))?;
            Ok((name, tensor))
        })
        .collect()
}

/// The entries of a zip archive whose members are all stored uncompressed, as `torch.save`
/// writes them (PyTorch's writer never compresses).
struct StoredZip {
    /// Member name to the offset of its local file header and its size.
    entries: HashMap<String, (usize, usize)>,
}

fn le_u16(b: &[u8], at: usize) -> Result<usize, String> {
    b.get(at..at + 2).map(|s| usize::from(u16::from_le_bytes([s[0], s[1]]))).ok_or_else(|| "truncated zip".into())
}

fn le_u32(b: &[u8], at: usize) -> Result<u64, String> {
    b.get(at..at + 4)
        .map(|s| u64::from(u32::from_le_bytes(s.try_into().expect("4 bytes"))))
        .ok_or_else(|| "truncated zip".into())
}

fn le_u64(b: &[u8], at: usize) -> Result<u64, String> {
    b.get(at..at + 8).map(|s| u64::from_le_bytes(s.try_into().expect("8 bytes"))).ok_or_else(|| "truncated zip".into())
}

fn to_usize(v: u64) -> Result<usize, String> {
    usize::try_from(v).map_err(|e| e.to_string())
}

impl StoredZip {
    const EOCD: &[u8] = b"PK\x05\x06";
    const ZIP64_LOCATOR: &[u8] = b"PK\x06\x07";
    const ZIP64_EOCD: &[u8] = b"PK\x06\x06";
    const CENTRAL: &[u8] = b"PK\x01\x02";
    const LOCAL: &[u8] = b"PK\x03\x04";

    fn parse(data: &[u8]) -> Result<Self, String> {
        // The end-of-central-directory record is the last 22 bytes plus an optional comment. A
        // file shorter than that (an empty or cut-off download) is no zip, not a panic.
        if data.len() < 22 {
            return Err(format!("not a zip archive ({} bytes; is the model file complete?)", data.len()));
        }
        let eocd = (0..=data.len() - 22)
            .rev()
            .take(22 + usize::from(u16::MAX))
            .find(|&i| {
                data.get(i..i + 4) == Some(Self::EOCD)
                    && i + 22 + le_u16(data, i + 20).unwrap_or(usize::MAX) == data.len()
            })
            .ok_or("not a zip archive (is the model file complete?)")?;
        let mut count = le_u16(data, eocd + 10)? as u64;
        let mut cd_offset = le_u32(data, eocd + 16)?;
        if let Some(loc) = eocd.checked_sub(20).filter(|&l| data.get(l..l + 4) == Some(Self::ZIP64_LOCATOR)) {
            let rec = to_usize(le_u64(data, loc + 8)?)?;
            if data.get(rec..rec + 4) != Some(Self::ZIP64_EOCD) {
                return Err("bad zip64 end record".into());
            }
            count = le_u64(data, rec + 32)?;
            cd_offset = le_u64(data, rec + 48)?;
        }

        let mut entries = HashMap::new();
        let mut p = to_usize(cd_offset)?;
        for _ in 0..count {
            if data.get(p..p + 4) != Some(Self::CENTRAL) {
                return Err(format!("bad central directory entry at byte {p}"));
            }
            let method = le_u16(data, p + 10)?;
            let mut size = le_u32(data, p + 20)?;
            let (name_len, extra_len, comment_len) =
                (le_u16(data, p + 28)?, le_u16(data, p + 30)?, le_u16(data, p + 32)?);
            let mut offset = le_u32(data, p + 42)?;
            let name = data.get(p + 46..p + 46 + name_len).ok_or("truncated zip")?;
            let name = String::from_utf8(name.to_vec()).map_err(|e| e.to_string())?;
            // A zip64 extra field (id 1) holds, in order, whichever of the uncompressed size,
            // compressed size and header offset overflowed 32 bits.
            let extra = data.get(p + 46 + name_len..p + 46 + name_len + extra_len).ok_or("truncated zip")?;
            let mut e = 0;
            while e + 4 <= extra.len() {
                let (id, len) = (le_u16(extra, e)?, le_u16(extra, e + 2)?);
                if id == 1 {
                    let mut f = e + 4;
                    for field in [&mut size, &mut 0u64, &mut offset] {
                        if *field == u64::from(u32::MAX) && f + 8 <= e + 4 + len {
                            *field = le_u64(extra, f)?;
                            f += 8;
                        }
                    }
                }
                e += 4 + len;
            }
            if method != 0 {
                return Err(format!("{name} is compressed (method {method})"));
            }
            entries.insert(name, (to_usize(offset)?, to_usize(size)?));
            p += 46 + name_len + extra_len + comment_len;
        }
        Ok(Self { entries })
    }

    fn get<'d>(&self, data: &'d [u8], name: &str) -> Result<&'d [u8], String> {
        let &(header, size) = self.entries.get(name).ok_or("missing from the archive")?;
        if data.get(header..header + 4) != Some(Self::LOCAL) {
            return Err("bad local file header".into());
        }
        let start = header + 30 + le_u16(data, header + 26)? + le_u16(data, header + 28)?;
        data.get(start..start + size).ok_or_else(|| "truncated member".into())
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum StorageType {
    Float,
    Long,
}

#[derive(Debug, Clone, PartialEq, Eq)]
struct StorageRef {
    key: String,
    kind: StorageType,
}

#[derive(Debug, Clone, PartialEq, Eq)]
struct TensorSpec {
    storage: StorageRef,
    offset: usize,
    shape: Vec<usize>,
    stride: Vec<usize>,
}

impl TensorSpec {
    /// Gathers the strided view out of the storage bytes into a contiguous tensor.
    fn materialize(&self, bytes: &[u8]) -> Result<CheckpointTensor, String> {
        let numel: usize = self.shape.iter().product();
        let mut indices = Vec::with_capacity(numel);
        let mut idx = vec![0usize; self.shape.len()];
        for _ in 0..numel {
            indices.push(self.offset + idx.iter().zip(&self.stride).map(|(i, s)| i * s).sum::<usize>());
            for d in (0..idx.len()).rev() {
                idx[d] += 1;
                if idx[d] < self.shape[d] {
                    break;
                }
                idx[d] = 0;
            }
        }
        let width = match self.storage.kind {
            StorageType::Float => 4,
            StorageType::Long => 8,
        };
        let elem =
            |i: usize| bytes.get(i * width..(i + 1) * width).ok_or_else(|| format!("element {i} out of storage"));
        let values = match self.storage.kind {
            StorageType::Float => TensorValues::F32(
                indices
                    .iter()
                    .map(|&i| Ok(f32::from_le_bytes(elem(i)?.try_into().expect("4 bytes"))))
                    .collect::<Result<_, String>>()?,
            ),
            StorageType::Long => TensorValues::I64(
                indices
                    .iter()
                    .map(|&i| Ok(i64::from_le_bytes(elem(i)?.try_into().expect("8 bytes"))))
                    .collect::<Result<_, String>>()?,
            ),
        };
        Ok(CheckpointTensor { shape: self.shape.clone(), values })
    }
}

/// Values on the unpickler stack.
#[derive(Debug, Clone)]
enum Value {
    Mark,
    None,
    /// A bool or float (`requires_grad`, metadata); its value is not needed.
    Scalar,
    Int(i64),
    Str(String),
    Tuple(Vec<Value>),
    List(Vec<Value>),
    Dict(Vec<(Value, Value)>),
    Global(String, String),
    Storage(StorageRef),
    Tensor(TensorSpec),
    /// Any other object; only its identity matters (e.g. backward-hook dicts).
    Opaque,
}

struct Unpickler<'a> {
    data: &'a [u8],
    pos: usize,
    stack: Vec<Value>,
    memo: HashMap<u32, Value>,
}

impl<'a> Unpickler<'a> {
    fn new(data: &'a [u8]) -> Self {
        Self { data, pos: 0, stack: Vec::new(), memo: HashMap::new() }
    }

    fn take(&mut self, n: usize) -> Result<&'a [u8], String> {
        let bytes = self.data.get(self.pos..self.pos + n).ok_or("truncated pickle")?;
        self.pos += n;
        Ok(bytes)
    }

    fn u8(&mut self) -> Result<u8, String> {
        Ok(self.take(1)?[0])
    }

    fn u32(&mut self) -> Result<u32, String> {
        Ok(u32::from_le_bytes(self.take(4)?.try_into().expect("4 bytes")))
    }

    fn line(&mut self) -> Result<String, String> {
        let rest = &self.data[self.pos..];
        let end = rest.iter().position(|&b| b == b'\n').ok_or("unterminated GLOBAL")?;
        self.pos += end + 1;
        Ok(String::from_utf8_lossy(&rest[..end]).into_owned())
    }

    fn pop(&mut self) -> Result<Value, String> {
        self.stack.pop().ok_or_else(|| "pickle stack underflow".to_owned())
    }

    fn pop_mark(&mut self) -> Result<Vec<Value>, String> {
        let mark = self.stack.iter().rposition(|v| matches!(v, Value::Mark)).ok_or("MARK not found")?;
        let items = self.stack.split_off(mark + 1);
        self.stack.pop();
        Ok(items)
    }

    fn top_dict(&mut self) -> Result<&mut Vec<(Value, Value)>, String> {
        match self.stack.last_mut() {
            Some(Value::Dict(items)) => Ok(items),
            Some(Value::Opaque) => Err("SETITEM on an unsupported object".into()),
            _ => Err("SETITEM without a dict".into()),
        }
    }

    /// Runs the pickle and returns the state dict's `(name, tensor)` entries.
    fn run(mut self) -> Result<Vec<(String, TensorSpec)>, String> {
        loop {
            let op = self.u8()?;
            match op {
                0x80 => {
                    let proto = self.u8()?;
                    if proto > 3 {
                        return Err(format!("pickle protocol {proto}"));
                    }
                }
                b'c' => {
                    let module = self.line()?;
                    let name = self.line()?;
                    self.stack.push(Value::Global(module, name));
                }
                b'q' => {
                    let i = u32::from(self.u8()?);
                    self.memo.insert(i, self.stack.last().ok_or("BINPUT on empty stack")?.clone());
                }
                b'r' => {
                    let i = self.u32()?;
                    self.memo.insert(i, self.stack.last().ok_or("LONG_BINPUT on empty stack")?.clone());
                }
                b'h' => {
                    let i = u32::from(self.u8()?);
                    self.stack.push(self.memo.get(&i).ok_or("BINGET of unknown memo")?.clone());
                }
                b'j' => {
                    let i = self.u32()?;
                    self.stack.push(self.memo.get(&i).ok_or("LONG_BINGET of unknown memo")?.clone());
                }
                b'(' => self.stack.push(Value::Mark),
                b')' => self.stack.push(Value::Tuple(Vec::new())),
                b't' => {
                    let items = self.pop_mark()?;
                    self.stack.push(Value::Tuple(items));
                }
                0x85..=0x87 => {
                    let n = usize::from(op - 0x84);
                    let at = self.stack.len().checked_sub(n).ok_or("TUPLEn underflow")?;
                    let items = self.stack.split_off(at);
                    self.stack.push(Value::Tuple(items));
                }
                b'}' => self.stack.push(Value::Dict(Vec::new())),
                b']' => self.stack.push(Value::List(Vec::new())),
                b'N' => self.stack.push(Value::None),
                0x88 | 0x89 => self.stack.push(Value::Scalar),
                b'K' => {
                    let v = self.u8()?;
                    self.stack.push(Value::Int(i64::from(v)));
                }
                b'M' => {
                    let v = u16::from_le_bytes(self.take(2)?.try_into().expect("2 bytes"));
                    self.stack.push(Value::Int(i64::from(v)));
                }
                b'J' => {
                    let v = i32::from_le_bytes(self.take(4)?.try_into().expect("4 bytes"));
                    self.stack.push(Value::Int(i64::from(v)));
                }
                0x8a => {
                    let n = usize::from(self.u8()?);
                    let bytes = self.take(n)?;
                    if n > 8 {
                        return Err("LONG1 wider than 64 bits".into());
                    }
                    let fill = if bytes.last().is_some_and(|b| b & 0x80 != 0) { 0xff } else { 0 };
                    let mut buf = [fill; 8];
                    buf[..n].copy_from_slice(bytes);
                    self.stack.push(Value::Int(i64::from_le_bytes(buf)));
                }
                b'G' => {
                    self.take(8)?;
                    self.stack.push(Value::Scalar);
                }
                b'X' => {
                    let n = self.u32()? as usize;
                    let s = String::from_utf8(self.take(n)?.to_vec()).map_err(|e| e.to_string())?;
                    self.stack.push(Value::Str(s));
                }
                b'a' => {
                    let item = self.pop()?;
                    match self.stack.last_mut() {
                        Some(Value::List(items)) => items.push(item),
                        _ => return Err("APPEND without a list".into()),
                    }
                }
                b'e' => {
                    let items = self.pop_mark()?;
                    match self.stack.last_mut() {
                        Some(Value::List(list)) => list.extend(items),
                        _ => return Err("APPENDS without a list".into()),
                    }
                }
                b's' => {
                    let value = self.pop()?;
                    let key = self.pop()?;
                    self.top_dict()?.push((key, value));
                }
                b'u' => {
                    let items = self.pop_mark()?;
                    if items.len() % 2 != 0 {
                        return Err("SETITEMS with an odd item count".into());
                    }
                    let mut it = items.into_iter();
                    let dict = self.top_dict()?;
                    while let (Some(k), Some(v)) = (it.next(), it.next()) {
                        dict.push((k, v));
                    }
                }
                b'Q' => {
                    let pid = self.pop()?;
                    self.stack.push(persistent_load(pid)?);
                }
                b'R' => {
                    let args = self.pop()?;
                    let callable = self.pop()?;
                    self.stack.push(reduce(callable, args)?);
                }
                b'b' => {
                    // BUILD sets object state, e.g. the state dict's `_metadata`; not needed.
                    self.pop()?;
                }
                b'.' => break,
                other => return Err(format!("unsupported pickle opcode 0x{other:02x} at byte {}", self.pos - 1)),
            }
        }
        let Value::Dict(mut items) = self.pop()? else { return Err("checkpoint is not a dict".into()) };
        // `_load_torch_checkpoint` also accepts {"state_dict": {...}, ...}.
        if let Some(pos) = items.iter().position(|(k, _)| matches!(k, Value::Str(s) if s == "state_dict")) {
            let Value::Dict(inner) = items.swap_remove(pos).1 else { return Err("state_dict is not a dict".into()) };
            items = inner;
        }
        items
            .into_iter()
            .map(|(k, v)| match (k, v) {
                (Value::Str(name), Value::Tensor(spec)) => Ok((name, spec)),
                (k, _) => Err(format!("entry {k:?} is not a tensor")),
            })
            .collect()
    }
}

/// `persistent_load` of `torch.serialization`: `('storage', <type>, key, location, numel)`.
fn persistent_load(pid: Value) -> Result<Value, String> {
    let Value::Tuple(items) = pid else { return Err("persistent id is not a tuple".into()) };
    match items.as_slice() {
        [Value::Str(tag), Value::Global(module, ty), Value::Str(key), _, _]
            if tag == "storage" && module == "torch" =>
        {
            let kind = match ty.as_str() {
                "FloatStorage" => StorageType::Float,
                "LongStorage" => StorageType::Long,
                other => return Err(format!("unsupported storage type {other}")),
            };
            Ok(Value::Storage(StorageRef { key: key.clone(), kind }))
        }
        _ => Err(format!("unsupported persistent id {items:?}")),
    }
}

fn usizes(v: &Value) -> Result<Vec<usize>, String> {
    let Value::Tuple(items) = v else { return Err("expected a tuple of ints".into()) };
    items
        .iter()
        .map(|i| match i {
            Value::Int(n) => usize::try_from(*n).map_err(|e| e.to_string()),
            _ => Err("expected an int".into()),
        })
        .collect()
}

fn reduce(callable: Value, args: Value) -> Result<Value, String> {
    let Value::Global(module, name) = &callable else { return Err("REDUCE of a non-global".into()) };
    match (module.as_str(), name.as_str()) {
        ("collections", "OrderedDict") => Ok(Value::Dict(Vec::new())),
        ("torch._utils", "_rebuild_tensor_v2") => {
            let Value::Tuple(a) = args else { return Err("tensor args are not a tuple".into()) };
            let [Value::Storage(storage), Value::Int(offset), shape, stride, ..] = a.as_slice() else {
                return Err("unexpected _rebuild_tensor_v2 arguments".into());
            };
            let (shape, stride) = (usizes(shape)?, usizes(stride)?);
            if shape.len() != stride.len() {
                return Err("shape and stride differ in length".into());
            }
            let offset = usize::try_from(*offset).map_err(|e| e.to_string())?;
            Ok(Value::Tensor(TensorSpec { storage: storage.clone(), offset, shape, stride }))
        }
        _ => Ok(Value::Opaque),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn spec(offset: usize, shape: Vec<usize>, stride: Vec<usize>) -> TensorSpec {
        TensorSpec { storage: StorageRef { key: "0".into(), kind: StorageType::Float }, offset, shape, stride }
    }

    fn floats(n: usize) -> Vec<u8> {
        (0..n).flat_map(|i| (i as f32).to_le_bytes()).collect()
    }

    #[test]
    fn contiguous_view() {
        let t = spec(0, vec![2, 3], vec![3, 1]).materialize(&floats(6)).unwrap();
        assert_eq!(t.values, TensorValues::F32(vec![0., 1., 2., 3., 4., 5.]));
    }

    #[test]
    fn transposed_view_with_offset() {
        // A 3x2 transpose of a 2x3 block starting at element 1.
        let t = spec(1, vec![3, 2], vec![1, 3]).materialize(&floats(8)).unwrap();
        assert_eq!(t.values, TensorValues::F32(vec![1., 4., 2., 5., 3., 6.]));
    }

    #[test]
    fn a_short_or_cut_off_file_is_an_error_not_a_panic() {
        // An empty or half-downloaded model panicked slicing past its end.
        for len in [0usize, 1, 3, 4, 21, 22, 100] {
            let mut data = b"PK\x03\x04".repeat(len.div_ceil(4));
            data.truncate(len);
            assert!(StoredZip::parse(&data).is_err(), "{len} bytes");
        }
        let dir = std::env::temp_dir().join(format!("snv-infer-empty-{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let path = dir.join("pileup.pt");
        std::fs::write(&path, "").unwrap();
        let err = read_state_dict(&path).expect_err("an empty model");
        assert!(err.to_string().contains("not a zip archive"), "{err}");
        let _ = std::fs::remove_dir_all(&dir);
    }

    #[test]
    fn scalar_and_out_of_range() {
        let t = spec(2, vec![], vec![]).materialize(&floats(3)).unwrap();
        assert_eq!(t.values, TensorValues::F32(vec![2.]));
        assert!(spec(0, vec![4], vec![1]).materialize(&floats(3)).is_err());
    }
}
