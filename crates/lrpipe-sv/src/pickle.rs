//! The subset of Python's pickle format (protocol 2-5) that Sniffles2 `.snf` blocks use.
//!
//! Decoding builds a generic [`Value`] tree; objects keep their class name and state, so SNF
//! files written by Sniffles2 can be read without Python. Encoding writes protocol 4 pickles
//! that Python (and Sniffles2) can load.

use std::fmt;
use std::rc::Rc;

/// A decoded Python value.
#[derive(Debug, Clone, PartialEq)]
pub enum Value {
    None,
    Bool(bool),
    Int(i64),
    Float(f64),
    /// Shared, as the unpickler memo refers to strings again.
    Str(Rc<str>),
    Bytes(Vec<u8>),
    List(Vec<Value>),
    Tuple(Vec<Value>),
    /// Key-value pairs in insertion order.
    Dict(Vec<(Value, Value)>),
    Set(Vec<Value>),
    /// A class reference (`module.name`).
    Global(String, String),
    /// An instance created with `NEWOBJ` or `REDUCE` and filled by `BUILD`.
    Object {
        class: Box<Value>,
        args: Vec<Value>,
        state: Box<Value>,
    },
}

impl Value {
    pub fn as_int(&self) -> Option<i64> {
        match self {
            Value::Int(i) => Some(*i),
            Value::Bool(b) => Some(i64::from(*b)),
            _ => None,
        }
    }

    pub fn as_float(&self) -> Option<f64> {
        match self {
            Value::Float(f) => Some(*f),
            Value::Int(i) => Some(*i as f64),
            _ => None,
        }
    }

    pub fn as_str(&self) -> Option<&str> {
        match self {
            Value::Str(s) => Some(s),
            _ => None,
        }
    }

    pub fn as_bool(&self) -> Option<bool> {
        match self {
            Value::Bool(b) => Some(*b),
            Value::Int(i) => Some(*i != 0),
            _ => None,
        }
    }

    /// Looks up a string key in a dict.
    pub fn get(&self, key: &str) -> Option<&Value> {
        match self {
            Value::Dict(items) => items.iter().find(|(k, _)| k.as_str() == Some(key)).map(|(_, v)| v),
            _ => None,
        }
    }

    /// An object's attribute (from its state dict).
    pub fn attr(&self, name: &str) -> Option<&Value> {
        match self {
            Value::Object { state, .. } => state.get(name),
            _ => None,
        }
    }

    pub fn is_none(&self) -> bool {
        matches!(self, Value::None)
    }

    /// Looks up a string key in a dict, for changing its value.
    pub fn get_mut(&mut self, key: &str) -> Option<&mut Value> {
        match self {
            Value::Dict(items) => items.iter_mut().find(|(k, _)| k.as_str() == Some(key)).map(|(_, v)| v),
            _ => None,
        }
    }

    /// Sets a string key's value in a dict that has it.
    pub fn set(&mut self, key: &str, value: Value) {
        if let Some(v) = self.get_mut(key) {
            *v = value;
        }
    }

    /// An object's state dict, for changing its attributes.
    pub fn state_mut(&mut self) -> Option<&mut Value> {
        match self {
            Value::Object { state, .. } => Some(state),
            _ => None,
        }
    }
}

#[derive(Debug)]
pub struct PickleError(String);

impl fmt::Display for PickleError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "invalid pickle: {}", self.0)
    }
}

impl std::error::Error for PickleError {}

fn err<T>(msg: impl Into<String>) -> Result<T, PickleError> {
    Err(PickleError(msg.into()))
}

/// Decodes a pickle.
///
/// # Errors
///
/// Fails on truncated data and on opcodes outside the supported subset.
pub fn decode(data: &[u8]) -> Result<Value, PickleError> {
    let mut stack: Vec<Value> = Vec::new();
    let mut marks: Vec<usize> = Vec::new();
    let mut memo: Vec<Value> = Vec::new();
    let mut p = 0usize;
    let take = |p: &mut usize, n: usize| -> Result<&[u8], PickleError> {
        let s = data.get(*p..*p + n).ok_or_else(|| PickleError("truncated".into()))?;
        *p += n;
        Ok(s)
    };
    let pop = |stack: &mut Vec<Value>| stack.pop().ok_or_else(|| PickleError("stack underflow".into()));
    let pop_mark = |stack: &mut Vec<Value>, marks: &mut Vec<usize>| -> Result<Vec<Value>, PickleError> {
        let m = marks.pop().ok_or_else(|| PickleError("missing MARK".into()))?;
        if m > stack.len() {
            return err("bad MARK");
        }
        Ok(stack.split_off(m))
    };
    loop {
        let op = *data.get(p).ok_or_else(|| PickleError("missing STOP".into()))?;
        p += 1;
        match op {
            0x80 => {
                take(&mut p, 1)?; // PROTO
            }
            0x95 => {
                take(&mut p, 8)?; // FRAME
            }
            b'.' => return pop(&mut stack),  // STOP
            b'(' => marks.push(stack.len()), // MARK
            b'N' => stack.push(Value::None),
            0x88 => stack.push(Value::Bool(true)),
            0x89 => stack.push(Value::Bool(false)),
            b'J' => stack.push(Value::Int(i64::from(i32::from_le_bytes(take(&mut p, 4)?.try_into().expect("4"))))),
            b'K' => stack.push(Value::Int(i64::from(take(&mut p, 1)?[0]))),
            b'M' => stack.push(Value::Int(i64::from(u16::from_le_bytes(take(&mut p, 2)?.try_into().expect("2"))))),
            0x8a => {
                // LONG1: little-endian two's complement
                let n = take(&mut p, 1)?[0] as usize;
                let b = take(&mut p, n)?;
                if n > 8 {
                    return err("integer too large");
                }
                let mut v: i64 = 0;
                for (i, &x) in b.iter().enumerate() {
                    v |= i64::from(x) << (8 * i);
                }
                if n > 0 && n < 8 && b[n - 1] & 0x80 != 0 {
                    v -= 1i64 << (8 * n);
                }
                stack.push(Value::Int(v));
            }
            b'G' => stack.push(Value::Float(f64::from_be_bytes(take(&mut p, 8)?.try_into().expect("8")))),
            0x8c => {
                let n = take(&mut p, 1)?[0] as usize;
                stack.push(Value::Str(String::from_utf8_lossy(take(&mut p, n)?).into()));
            }
            b'X' => {
                let n = u32::from_le_bytes(take(&mut p, 4)?.try_into().expect("4")) as usize;
                stack.push(Value::Str(String::from_utf8_lossy(take(&mut p, n)?).into()));
            }
            0x8d => {
                let n = u64::from_le_bytes(take(&mut p, 8)?.try_into().expect("8")) as usize;
                stack.push(Value::Str(String::from_utf8_lossy(take(&mut p, n)?).into()));
            }
            b'C' => {
                let n = take(&mut p, 1)?[0] as usize;
                stack.push(Value::Bytes(take(&mut p, n)?.to_vec()));
            }
            b'B' => {
                let n = u32::from_le_bytes(take(&mut p, 4)?.try_into().expect("4")) as usize;
                stack.push(Value::Bytes(take(&mut p, n)?.to_vec()));
            }
            b']' => stack.push(Value::List(Vec::new())),
            b'}' => stack.push(Value::Dict(Vec::new())),
            b')' => stack.push(Value::Tuple(Vec::new())),
            0x8f => stack.push(Value::Set(Vec::new())),
            0x85 => {
                let a = pop(&mut stack)?;
                stack.push(Value::Tuple(vec![a]));
            }
            0x86 => {
                let b = pop(&mut stack)?;
                let a = pop(&mut stack)?;
                stack.push(Value::Tuple(vec![a, b]));
            }
            0x87 => {
                let c = pop(&mut stack)?;
                let b = pop(&mut stack)?;
                let a = pop(&mut stack)?;
                stack.push(Value::Tuple(vec![a, b, c]));
            }
            b't' => {
                let items = pop_mark(&mut stack, &mut marks)?;
                stack.push(Value::Tuple(items));
            }
            b'l' => {
                let items = pop_mark(&mut stack, &mut marks)?;
                stack.push(Value::List(items));
            }
            b'a' => {
                let v = pop(&mut stack)?;
                match stack.last_mut() {
                    Some(Value::List(l)) => l.push(v),
                    _ => return err("APPEND to non-list"),
                }
            }
            b'e' => {
                let items = pop_mark(&mut stack, &mut marks)?;
                match stack.last_mut() {
                    Some(Value::List(l)) => l.extend(items),
                    _ => return err("APPENDS to non-list"),
                }
            }
            b's' => {
                let v = pop(&mut stack)?;
                let k = pop(&mut stack)?;
                match stack.last_mut() {
                    Some(Value::Dict(d)) => set_item(d, k, v),
                    _ => return err("SETITEM on non-dict"),
                }
            }
            b'u' => {
                let items = pop_mark(&mut stack, &mut marks)?;
                match stack.last_mut() {
                    Some(Value::Dict(d)) => {
                        let mut it = items.into_iter();
                        while let (Some(k), Some(v)) = (it.next(), it.next()) {
                            set_item(d, k, v);
                        }
                    }
                    _ => return err("SETITEMS on non-dict"),
                }
            }
            0x90 => {
                let items = pop_mark(&mut stack, &mut marks)?;
                match stack.last_mut() {
                    Some(Value::Set(s)) => s.extend(items),
                    _ => return err("ADDITEMS on non-set"),
                }
            }
            0x91 => {
                let items = pop_mark(&mut stack, &mut marks)?;
                stack.push(Value::Set(items));
            }
            0x93 => {
                let name = pop(&mut stack)?;
                let module = pop(&mut stack)?;
                match (module, name) {
                    (Value::Str(m), Value::Str(n)) => stack.push(Value::Global(m.to_string(), n.to_string())),
                    _ => return err("STACK_GLOBAL needs two strings"),
                }
            }
            b'c' => {
                let line = |p: &mut usize| -> Result<String, PickleError> {
                    let rest = &data[*p..];
                    let end = rest.iter().position(|&b| b == b'\n').ok_or_else(|| PickleError("GLOBAL".into()))?;
                    *p += end + 1;
                    Ok(String::from_utf8_lossy(&rest[..end]).into_owned())
                };
                let m = line(&mut p)?;
                let n = line(&mut p)?;
                stack.push(Value::Global(m, n));
            }
            0x81 | b'R' => {
                // NEWOBJ (cls, args) / REDUCE (callable, args)
                let args = pop(&mut stack)?;
                let class = pop(&mut stack)?;
                let args = match args {
                    Value::Tuple(a) => a,
                    _ => return err("object arguments must be a tuple"),
                };
                stack.push(Value::Object { class: Box::new(class), args, state: Box::new(Value::None) });
            }
            b'b' => {
                let st = pop(&mut stack)?;
                match stack.last_mut() {
                    Some(Value::Object { state, .. }) => **state = st,
                    _ => return err("BUILD on non-object"),
                }
            }
            0x94 => {
                // MEMOIZE
                let v = stack.last().ok_or_else(|| PickleError("MEMOIZE on empty stack".into()))?.clone();
                memo.push(v);
            }
            b'q' => {
                let id = u32::from(take(&mut p, 1)?[0]);
                memoize(&mut memo, id, stack.last().cloned().ok_or_else(|| PickleError("BINPUT".into()))?);
            }
            b'r' => {
                let id = u32::from_le_bytes(take(&mut p, 4)?.try_into().expect("4"));
                memoize(&mut memo, id, stack.last().cloned().ok_or_else(|| PickleError("LONG_BINPUT".into()))?);
            }
            b'h' => {
                let id = u32::from(take(&mut p, 1)?[0]);
                stack.push(memo.get(id as usize).cloned().ok_or_else(|| PickleError("BINGET".into()))?);
            }
            b'j' => {
                let id = u32::from_le_bytes(take(&mut p, 4)?.try_into().expect("4"));
                stack.push(memo.get(id as usize).cloned().ok_or_else(|| PickleError("LONG_BINGET".into()))?);
            }
            0x30 => {
                pop(&mut stack)?; // POP
            }
            other => return err(format!("unsupported opcode 0x{other:02x} at byte {}", p - 1)),
        }
    }
}

/// `dict[k] = v`. Python writes each key of a dict once, so keys are appended without a search.
fn set_item(d: &mut Vec<(Value, Value)>, k: Value, v: Value) {
    d.push((k, v));
}

/// Stores a value under an explicit memo index.
fn memoize(memo: &mut Vec<Value>, id: u32, v: Value) {
    let id = id as usize;
    if memo.len() <= id {
        memo.resize(id + 1, Value::None);
    }
    memo[id] = v;
}

/// Encodes a value as a protocol 4 pickle.
pub fn encode(v: &Value) -> Vec<u8> {
    let mut out = vec![0x80, 4];
    write(&mut out, v);
    out.push(b'.');
    out
}

fn write(out: &mut Vec<u8>, v: &Value) {
    match v {
        Value::None => out.push(b'N'),
        Value::Bool(true) => out.push(0x88),
        Value::Bool(false) => out.push(0x89),
        Value::Int(i) => {
            if (0..=0xFF).contains(i) {
                out.push(b'K');
                out.push(*i as u8);
            } else if (0..=0xFFFF).contains(i) {
                out.push(b'M');
                out.extend_from_slice(&(*i as u16).to_le_bytes());
            } else if let Ok(i32v) = i32::try_from(*i) {
                out.push(b'J');
                out.extend_from_slice(&i32v.to_le_bytes());
            } else {
                out.push(0x8a);
                out.push(8);
                out.extend_from_slice(&i.to_le_bytes());
            }
        }
        Value::Float(f) => {
            out.push(b'G');
            out.extend_from_slice(&f.to_be_bytes());
        }
        Value::Str(s) => {
            if s.len() < 256 {
                out.push(0x8c);
                out.push(s.len() as u8);
            } else {
                out.push(b'X');
                out.extend_from_slice(&(s.len() as u32).to_le_bytes());
            }
            out.extend_from_slice(s.as_bytes());
        }
        Value::Bytes(b) => {
            out.push(b'B');
            out.extend_from_slice(&(b.len() as u32).to_le_bytes());
            out.extend_from_slice(b);
        }
        Value::List(items) => {
            out.push(b']');
            if !items.is_empty() {
                out.push(b'(');
                items.iter().for_each(|x| write(out, x));
                out.push(b'e');
            }
        }
        Value::Tuple(items) => match items.len() {
            0 => out.push(b')'),
            1..=3 => {
                items.iter().for_each(|x| write(out, x));
                out.push([0x85, 0x86, 0x87][items.len() - 1]);
            }
            _ => {
                out.push(b'(');
                items.iter().for_each(|x| write(out, x));
                out.push(b't');
            }
        },
        Value::Dict(items) => {
            out.push(b'}');
            if !items.is_empty() {
                out.push(b'(');
                for (k, x) in items {
                    write(out, k);
                    write(out, x);
                }
                out.push(b'u');
            }
        }
        Value::Set(items) => {
            out.push(0x8f);
            if !items.is_empty() {
                out.push(b'(');
                items.iter().for_each(|x| write(out, x));
                out.push(0x90);
            }
        }
        Value::Global(m, n) => {
            write(out, &Value::Str(m.as_str().into()));
            write(out, &Value::Str(n.as_str().into()));
            out.push(0x93);
        }
        Value::Object { class, args, state } => {
            write(out, class);
            write(out, &Value::Tuple(args.clone()));
            out.push(0x81);
            if !state.is_none() {
                write(out, state);
                out.push(b'b');
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn round_trips() {
        let v = Value::Dict(vec![
            (Value::Str("INS".into()), Value::List(vec![Value::Int(1), Value::Int(-70000), Value::Int(1 << 40)])),
            (Value::Int(600000), Value::Tuple(vec![Value::Int(0), Value::Str(".".into()), Value::None])),
            (Value::Str("f".into()), Value::Float(11.523888232710346)),
            (
                Value::Str("obj".into()),
                Value::Object {
                    class: Box::new(Value::Global("sniffles.sv".into(), "SVCall".into())),
                    args: vec![],
                    state: Box::new(Value::Dict(vec![(Value::Str("pos".into()), Value::Int(641939))])),
                },
            ),
            (Value::Str("b".into()), Value::Bool(true)),
        ]);
        assert_eq!(decode(&encode(&v)).expect("valid pickle"), v);
    }

    #[test]
    fn decodes_negative_long1() {
        // LONG1 with 2 bytes 0x00 0xFF = -256
        assert_eq!(decode(&[0x80, 4, 0x8a, 2, 0x00, 0xFF, b'.']).expect("valid"), Value::Int(-256));
    }
}
