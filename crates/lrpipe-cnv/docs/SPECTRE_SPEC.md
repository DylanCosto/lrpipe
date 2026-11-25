# Spectre 0.2.1 (commit 90580943): implementation-level specification of `spectre CNVCaller`

Source tree: Spectre 0.2.1, commit 90580943 (https://github.com/fritzsedlazeck/Spectre) (paths below are relative to it). Every rule cites `file:line`.
"py" = Python built-in semantics, "np" = numpy, "pd" = pandas.

## 0. Reference environment (what "exact parity" means)

* Python 3.10 / 3.11 (README.md:17, Dockerfile:1 `python:3.11.4`, CI `.github/workflows/release.yaml:13` 3.10).
* numpy **1.x** only (1.24.3 pinned in README.md:34). numpy 2.x breaks Spectre: `np.NaN` (coverage_stats.py:11-18,37-41; analysis.py:243) and `np.msort` (unidip.py:38, dip.py:76) were removed in 2.0, so a numpy-2 run crashes when the first `CoverageStatistics` is built.
* pandas **2.0.x** (README.md:34 pins 2.0.1). Some float behaviour depends on the pandas major version; see §5.4.
* scipy 1.10.1 (`scipy.stats.norm.sf`, `ks_2samp`, of which only `norm.sf` is reached).
* pysam 0.22.0, which bundles htslib/bcftools 1.18. The output VCF is re-serialized by bcftools (§4.2).
* Bottleneck **not** installed. If it is installed, pandas reductions change summation order (§5.4).
* Version strings: logged `0.2.1-patch-july15` (spectre.py:157); written into `.spc` as `"0.2.alpha"` (analysis.py:810).

---

## 1. CLI

Entry point `spectre = spectre.spectre:main` (setup.cfg:22-24). Parser: spectre.py:405-630. The `main()` dispatcher is at spectre.py:633-673.

### 1.1 Subcommands
| subcommand | lines | action |
|---|---|---|
| `CNVCaller` | spectre.py:496-574, 652-655 | single-sample calling (this spec) |
| `RemoveNs` | spectre.py:582-598, 656-660 | build `.mdr` from FASTA (§2.3) |
| `Population` | spectre.py:603-621, 662-665 | population mode (§8) |
| `version` | spectre.py:490-491, 666-667 | logs `Spectre version: 0.2.1-patch-july15` to stderr |
| `dev` | spectre.py:626 | no options; falls through to `else` and prints help (spectre.py:668-670) |

argparse errors exit with code 2. When no subcommand is given, `len(sys.argv) < 2` prints help and exits 1 (spectre.py:641-643); this check runs *after* `parse_args`. An unknown or empty command prints help. The program always ends with `sys.exit(0)` (spectre.py:673) unless an earlier `sys.exit(1)` or an exception fires. All logging goes to **stderr** with the format `lrpipe_cnv::%(asctime)s::%(levelname)s::%(name)s>  %(message)s` (logger.py:8). The level is DEBUG with `--dev`, otherwise INFO (logger.py:6).

### 1.2 `CNVCaller` options (spectre.py:500-574)
| flag(s) | dest | type | default | req | used for |
|---|---|---|---|---|---|
| `-c/--coverage` | coverage_dir | str, `nargs='+'` | — | **yes** | mosdepth `*.regions.bed.gz`; single mode uses `[0]` (spectre.py:324) |
| `-s/--sample-id` | sample_id | str, `nargs='+'` | — | **yes** | output names, RNG seed |
| `-d/--output-dir` | output_dir | str | `.` | **yes** | output dir (created, spectre.py:303-304) |
| `-r/--reference` | reference | str | — | **yes** | FASTA (needs `.fai`; pysam creates it if absent) |
| `-m/--metadata` | metadata_filename | str | `""` | no | `.mdr` path; created if missing (§2.3) |
| `-v/--snv` | snv_file | str | `""` | no | SNV VCF (§2.4) |
| `-sj/--snfj` | snfj_file | str | `""` | no | Sniffles SNFJ breakpoints (§3.12) |
| `-l/--blacklist` | black_list_file | str | `""` | no | blacklist BED (§2.3) |
| `-o/--only-chr` | only_chr_list | str | `""` | no | comma list; default = all genome chromosomes (spectre.py:280-281) |
| `-p/--ploidy` | ploidy | int | 2 | no | global ploidy |
| `-pc/--ploidy-chr` | ploidy_chr | str | `""` | no | `chr:N,...` (§2.5) |
| `-n/--n-size` | n_size | int | 5 | no | N-run length for `.mdr` creation |
| `-mcl/--min-cnv-len` | min_cnv_len | int | 100000 | no | final size filter **and** a startup check (§7 bug B1) |
| `-mchl/--min-chr-len` | min_chr_len | int | 1000000 | no | keep a FASTA contig iff `len > min_chr_len` (spectre.py:185) |
| `-t/--threads` | threads | int | 1 | no | population mode only |
| `-i/--population` | run_population_mode | flag | False | no | population mode |
| `-sco/--sample-coverage-overwrite` | sample_coverage_overwrite | float | None | no | replaces the genome median (analysis.py:168-170) |
| `-dmc/--disable-max-coverage` | disable_max_coverage | flag | False | no | skip outlier removal (analysis.py:239) |
| `-a/--cancer` | is_cancer | flag | False | no | disables threshold clamping (cnv_metrics.py:439) |
| `-qcna/--qc-noise-allow` | qc_noise_allow | float | 0.5 | no | clean_by_threshold (analysis.py:644-645) |
| `-lohkb/--loh-min-snv-perkb` | loh_min_snv_perkb | int | **10** | no | LoH filter (README says 5) |
| `-lohsnv/--loh-min-snv-total` | loh_min_snv_total | int | **1000** | no | LoH filter (README says 100) |
| `-lohsize/--loh-min-region-size` | loh_min_region_size | float | **1000000.0** | no | LoH filter (README says 100000) |
| `-lohonly/--loh-only` | loh_only | flag | False | no | skip CNV calling (crashes, §7 B2) |
| `-cn/--cn-neutral` | cn_neutral | float | 0.0 | no | only with `--loh-only` (spectreCNV.py:66-69) |
| `-0/--dev` | as_dev | flag | False | no | debug outputs + DEBUG log |
| `-01/--dev-max-std-outlier-rm` | max_std_outlier_rm | int | 5 | no | outlier cutoff (analysis.py:75, 242) |
| `-02/--mosdepth-cov-genome-chr-diff` | … | float | 0.10 | no | only written to `.spc` (analysis.py:821) |
| `-03/--lower-2n-threshold` | … | float | 1.5 | no | overwritten before use (analysis.py:211) — dead |
| `-04/--upper-2n-threshold` | … | float | 2.5 | no | overwritten (analysis.py:212) — dead |
| `-05/--cov-diff-threshold` | cov_diff_threshold | float | 0.80 | no | merge criterion (analysis.py:404, 416-417) |
| `-06/--dist-proportion` | dist_proportion | float | 0.25 | no | read but unused (analysis.py:406); written to `.spc` |

Hidden or fixed parameters: `min_chr_length=1e5`, `dist_min_overwrite=10000` (analysis.py:73-74), `min_run=10` (call_cnv_coverage.py:11), `bin_size` inferred (§2.1). The dev params end up as `spectre_args.dev_*` (spectre.py:87-92). Note that `SpectreCallParam.max_std_outlier_rm` etc. (spectre.py:42-48) keep their defaults and are **not** what the analysis reads.

### 1.3 Validation that actually happens
* argparse enforces the required flags and the int/float types.
* `--ploidy-chr`: each comma item must contain `:`, else it logs an error and calls `sys.exit(1)` (spectre.py:288-292). `chr, p = item.split(":")` raises ValueError on 2+ colons. A chromosome not among the genome chromosomes is skipped with a warning (spectre.py:295-297). `int(p)` raises ValueError on non-integers, and p=0 raises ZeroDivisionError.
* `--only-chr` is **not** validated. A name with no coverage or outside the genome raises KeyError at analysis.py:187.
* Multiple `--sample-id` without `--population` logs a warning and analyses only `sample_id[0]` with `coverage_dir[0]` (spectre.py:352-355, 324). Pairs are formed with `zip`, so the shorter list wins.
* `min_cnv_len` check, analysis.py:86-101, 146-147: the run continues **only if `10000 < min_cnv_len <= 100000`**; otherwise it calls `sys.exit(1)` (bug B1).
* SNV VCF: the header must declare FORMAT `GT`, `DP` and (`AF` or `VAF`), else `sys.exit(1)` (snv_analysis.py:71-75).
* Coverage: if the path ends with `.regions.bed.gz` and `<path>.csi` is missing, a CSI index is built (spectreCNV.py:40-48). Any other path with no `.csi` fails at `TabixFile` open (analysis.py:112).
* No coverage for any genome chromosome: error, then `sys.exit(1)` (analysis.py:155-159).

### 1.4 Path handling
`out_dir` and `reference` are made absolute (spectre.py:269-270). The coverage path is used as given (`coverage_dir[0]`, spectre.py:324; the abspath at line 323 is discarded). The SNV path is made absolute in `vcf_has_index` (vcf_parser.py:30). Directories created: `out_dir`, `out_dir/tmp`, `out_dir/img`, and `out_dir/debug` with `--dev` (spectre.py:303-318). `debug_dir` is `out_dir` itself when not in dev mode (spectre.py:318). An empty sample id becomes `"sample-"` (spectreCNV.py:16).

---

## 2. Inputs and parsing

### 2.1 mosdepth coverage
* Expected invocation (README.md:63): `mosdepth -t 8 -x -b 1000 -Q 20 <prefix> <bam>`, i.e. fixed 1 kb windows, MAPQ ≥ 20, fast mode. Spectre reads **only** `<prefix>.regions.bed.gz` plus its `.csi` index (spectre.py:411-413, analysis.py:112). It never reads the mosdepth summary or per-base files (`MosdepthReader` is dead code, mosdepthReader.py:26). `-n` (no per-base output) is harmless.
* Reading (analysis.py:117-131): for each `chrom` in `genome_info["chromosomes"]` (FASTA order, §2.2), skip it with a warning if it is not in the tabix index's contigs (analysis.py:119-121). Otherwise fetch all lines of that contig in file order and split each on `\t` into **exactly 4 fields** (`chrom, start, end, depth`; a 5-column file raises ValueError). Store `positions = int(start)` (0-based bin start) and `coverage = float(depth)` (python float parse; mosdepth prints `%.2f`). The `end` column is ignored.
* **Bin size** is `abs(positions[1] - positions[0])` of the *first chromosome that has coverage* (np.int64), computed once (analysis.py:137-138). This needs ≥2 lines on that contig (else IndexError). The metadata is loaded at the same moment (analysis.py:141-144).

### 2.2 Reference / genome info
`make_genome_info` (spectre.py:181-193) reads `pysam.FastaFile(reference)`. The reference is used for nothing except contig names and lengths, plus `.mdr` creation. For each `(name, len)` in `.fai` order, the contig is kept iff `len > min_chr_len` (strict). The result is `genome_info = {"chromosomes":[...], "chr_lengths":{name:len}}`, extended with `"chr_ploidy"` (spectre.py:284-300).

**Contig names** are never normalized: no `chr` prefix handling and no aliases. Names must match exactly between FASTA, mosdepth, `.mdr`, blacklist, SNV VCF and SNFJ.

### 2.3 Metadata: `.mdr` (N regions) + blacklist BED
**Which `.mdr` path** (spectre.py:196-225, called from 277). The path is `--metadata` if given, else `<out_dir>/metadata.mdr`, made absolute. If the file does not exist it is created with `FastaRef.get_n_regions` (§2.3.1) at exactly that path (metadataCollector.py:238-243; `os.path.join(out_dir, abs)` returns `abs`). If it exists it is used as-is.

#### 2.3.1 `.mdr` creation (`RemoveNs` and CNVCaller fallback), metadataCollector.py:25-99, 148-169
* FASTA open: `gzip.open(...,"rt")` if the string `"gz"` occurs anywhere in the path, else plain open (metadataCollector.py:46).
* For each line with `\n` stripped:
  * If the line starts with `>` **or `;`**, it is a header: `chromosome = line[1:].split(" ")[0]` (split on **space only**, not tabs). `ref_dict[chrom] = []` and `i = j = 0` (metadataCollector.py:79-89).
  * Otherwise, for each char `c` (upper-cased):
    * if `c == 'N'`: `j += 1`;
    * else: if `(j - i) >= n_size`, append `(i+1, j-1)` to `ref_dict[chrom]`; then `i = j; j += 1`. The non-N char also bumps the A/C/G/T counters (metadataCollector.py:51-77).
  * Semantics (quirk Q-MDR): with `j` = 0-based index of the current char and `i` = index of the previous non-N char (or 0 at contig start), a run of N at indices `i+1..j-1` is emitted when `j - i >= n_size`, i.e. **run length ≥ n_size − 1**. A leading N run starting at index 0 is reported as starting at 1. A trailing N run at the contig end is never reported. Coordinates are effectively 0-based inclusive `(first_N, last_N)`.
* Report text, written with `writelines` of one string (metadataCollector.py:150-166):
  ```
  #####\tStatistic report\n
  #####\tBasepair statistic\n
  BPDEF\tA\tT\tC\tG\tGC%\tTotal bp\n
  BPSTA\t{A}\t{T}\t{C}\t{G}\t{str((G+C)/j_last*100)}\t{j_last}\n
  #####\tN-Sequence positions in given reference file\n
  NSDEF\tFrom\tTo\n
  NSPOS\t{chrom}\t{start}\t{end}\n   (dict order = FASTA order, then list order)
  ```
  `j_last` is the `j` counter of the **last** contig only (bug: GC% uses the last contig's length; see data/grch38.mdr:4, which shows `16569` = chrM).

#### 2.3.2 Loading for analysis, `FastaRef.load_metadata(mdr, blacklist, bin_size)` (metadataCollector.py:129-146)
1. `.mdr` (metadataCollector.py:101-127): for each line starting with `NSPOS`, split on `\t`: `chrom=f[1]`, `start=int(f[2])`, `end=int(f[3])` (int() tolerates the trailing `\n`). Snap to the bin grid: `if start % b: start -= start % b`; `if end % b: end += b - end % b`. Append the **int** tuple to `nlist[chrom]` (dict insertion order = file order).
2. Blacklist (metadataCollector.py:194-214), only if a path is given. For each line, `term = line.decode()[:-1]` **removes the last char unconditionally**: on a final line without `\n` this chops a real character (data/grch37_blacklist_normal.bed and data/grch37_blacklist_plus.bed lack a trailing newline). Then `chrom, start, end = term.split("\t")[:3]`; the tuple is stored **as strings**, unsnapped. An empty line raises ValueError. Header or track lines are not skipped (they crash later in `int()`).
3. Merge (metadataCollector.py:217-236): `result = mdr.copy()` (shallow). For each blacklist chrom: if present, append every blacklist tuple **not already in the list** (string tuples never equal int tuples, so this only dedups exact duplicate blacklist lines); else `result[chrom] = blacklist_list`. Order: MDR tuples first, then blacklist tuples. Key order: MDR chroms, then new blacklist chroms.

The result `metadata: dict[chrom -> list[(start,end)]]` (int or str elements) is used in three places:
* **Masking** (analysis.py:293-308): for each tuple, `s=int(start)`, `e=int(end)`, and `coverage[(pos >= s) & (pos <= e)] = NaN`, where `pos` is the bin start. BED half-open semantics are ignored: blacklist `chr1 1 1000000` masks bins 1000…1000000 inclusive but **not** bin 0.
* **Merge blocking** (analysis.py:433-439), §3.6.
* **CNVMetrics exclusion zones** (cnv_metrics.py:65-80, 99-116), effectively a no-op (§3.4, Q-EXCL).

Shipped data: `data/grch38.mdr`, `data/grch37.mdr` and several blacklists (4 columns; one 3-column line at data/grch37_blacklist_plus.bed:134).

### 2.4 SNV VCF (`--snv`)
* Indexing (vcf_parser.py:29-44): if `<abs>.tbi` is missing, call `pysam.tabix_index(path, preset="vcf", force=True, keep_original=True)`. A plain `.vcf` is bgzipped to `.vcf.gz` (written next to the input), and the `.gz` path is then used. A file with only a `.csi` gets re-indexed. Then `pysam.VariantFile(path)` (snv_analysis.py:33).
* Tags (snv_analysis.py:67-79): `GT`, `DP`, and `af_tag = "AF"` if FORMAT `AF` is declared, else `"VAF"`. `gt_het=(0,1)`, `gt_ref=(0,0)`, `gt_alt=(1,1)` (pysam allele-index tuples; phased `0|1` is also `(0,1)`; `(1,0)`, `(1,2)` and `(None,None)` are neither het nor alt).
* **Value typing (pysam):** a FORMAT value is a scalar iff the header declares `Number=1` (fixed). For any other Number (`A`, `R`, `G`, `.`, `2`) it is a **tuple**, even with one ALT. Every reader skips records whose AF value is a tuple (snv_analysis.py:166-167, 217-218, 286-287; call_cnv_AF.py:90-92). So DeepVariant's `VAF` (Number=A) gets every record skipped, which leads to a crash at `np.max([])` (snv_analysis.py:104). Clair3's `AF` (Number=1) works.
* **Numeric precision:** AF and QUAL are stored as float32 by htslib and widened to float64 by pysam (e.g. `0.95` becomes `0.949999988079071`). DP is an int. A Rust port must parse AF and QUAL as f32 and widen to f64.
* Missing values: QUAL `.` gives `None`, and `None > 5` raises TypeError (crash). AF/DP `.` gives `None`, which crashes in `abs()`/`float()`/`int()`. Multi-sample: `[snv] = samples.values()` raises ValueError; the fallback at snv_analysis.py:162-163 is itself broken. Only single-sample VCFs work.
* Thresholds (snv_analysis.py:81-87): `max_deviation_perfect_het_af=0.12`, `perfect_het_af=0.5`, `min_snv_qual=5` (strict `qual > 5`), `min_coverage_site=3` (strict `DP > 3`, only in `snv_nhet_chromosome`), `min_hom_alt_af=0.75`, `min_het_chr_sites=0.4`, `min_gt_alt_prop=0.85`. `max_abs_dev_2x=0.2` (snv_analysis.py:65) is unused. LoH argument floors: perkb ≥ 1, total ≥ 10, region ≥ 10000.0 (snv_analysis.py:50-58).
* Region fetches use strings `f"{chrom}:{start}-{end}"`, which htslib reads as **1-based inclusive**. They return records whose `[POS-1, POS-1+rlen)` overlaps `[start-1, end)`, where `rlen = len(REF)` or derives from `INFO/END`. Indels whose REF spans into the region are therefore returned too.

### 2.5 Ploidy and sex chromosomes
`chr_ploidy[c] = ploidy` for every genome chromosome (spectre.py:284). For each `--ploidy-chr c:p`: `chr_ploidy[c] = int((1/(int(p)/ploidy))*ploidy)` = `int(ploidy²/p)` evaluated in float (spectre.py:299-300). With ploidy 2: `c:1` gives **4**, `c:2` gives 2, `c:3` gives 1, `c:4` gives 1. This value is a **coverage multiplier**: `normalized_cov_ploidy = coverage/norm * chr_ploidy[c]` (analysis.py:272-273). A haploid chrX is therefore scaled up to look diploid. Thresholds and clamps use the global `--ploidy` (analysis.py:204, cnv_metrics.py:441-449). There is no automatic sex detection and no chrX/chrY special-casing (data/grch38_y.txt is unused).

---

## 3. Algorithm (`SpectreCNV.cnv_call`, spectreCNV.py:59-134)

Order of steps:
1. if `--snv`: SNV CN-neutral estimation (§3.1);
2. unless `--loh-only`: `data_normalization` (§3.2–3.4), `call_cnv_coverage` (§3.5), `get_cnv_metrics()` (§3.10, unrefined), `refine_cnv_calls` (§3.6–3.9), `get_cnv_metrics(refined_cnvs=True)` (§3.10);
3. if `--snv`: AF concordance filter (§3.11, output-dead) and LoH (§3.13);
4. if `--snfj`: breakpoint support (§3.12);
5. CNV/LoH merge (§3.14), then `.spc`, BED, VCF and plots (§4).

### 3.1 SNV-derived CN-neutral coverage (snv_analysis.py:195-237, 151-193, 102-146)
Skipped (and `med=avg=--cn-neutral, std=1`) only when `cn_neutral != 0.0 and --loh-only` (spectreCNV.py:66-69).
1. `snv_copy_number_state(ploidy=2)` **forces `self.ploidy = 2`**, ignoring `--ploidy` (snv_analysis.py:200). This value is used later for LoH CN.
2. `snv_nhet_chromosome` (snv_analysis.py:151-193): iterate all records (`fetch()`, file order). `chr_all += 1`. Skip records whose AF is a tuple. If `qual > 5 and DP > 3`:
   * if GT == (0,1): when `abs(0.5 - AF) < 0.12`, do `chr_het += 1` and `cov_het += DP`;
   * else (any non-(0,1) GT): when `AF >= 0.75`, do `chr_alt += 1`.

   Keep chromosome `c` (iterating `genome_info["chromosomes"]`) iff `chr_alt[c] + chr_het[c] > 0` and `chr_het[c] / (chr_alt[c] + chr_het[c]) > 0.4`.
3. Second full pass (snv_analysis.py:210-226): for records on kept chromosomes, skip tuple AF; if `qual > 5`, GT == (0,1) and `abs(0.5 - AF) < 0.12`, append `DP` to `perfect_het_depth`. **There is no DP > 3 filter in this pass.** Order is file order.
4. `snv_multimodal_detect` (snv_analysis.py:102-128):
   * if `np.max(depths) > 10*np.nanmean(depths)`, apply `het_clean`: keep values `<= np.percentile(depths, 98)` (linear method; §6), preserving order. An empty list makes `np.max` raise ValueError (crash).
   * `UniDip(depths).run()` returns a list of index intervals `(s, e)` (§3.1.1).
   * With ≥1 interval: `(s, e) = intervals[0]`; `med = np.nanmedian(depths[s:e])`; `avg = np.nanmean(...)`; `std = np.nanstd(...)`. With 0 intervals, use the whole list.
   * **Quirk Q-UNIDIP:** the intervals index the *sorted* data (and, worse, the dip test's unique-value axis, §3.1.1), but they are applied to the **unsorted** VCF-order list, using `[s:e]` against an inclusive `e`. In the common unimodal case, `run()` returns `[(0, n-1)]`, so **`med = median(depths[0 : n-1])` = median of all het DPs except the last one in VCF order**.
5. Only `med` affects calls: it becomes the normalization value (analysis.py:231-234). `avg` and `std` are only logged.

#### 3.1.1 UniDip / dip test (unidip.py, dip.py), vendored from BenjaminDoran/unidip (GPL-3)
* `UniDip(dat, is_hist=False, alpha=0.05, ntrials=100, mrg_dst=1)`. `self.dat = np.msort(np.array(dat))` is a sorted copy (unidip.py:37-43).
* `run()` calls `merge_intervals(_unidip(0, len(dat), True))` (unidip.py:45-52). `merge_intervals` sorts the tuples and merges `idx[0] - last[1] <= 1` into `(last[0], idx[1])` (unidip.py:54-68).
* `_unidip(start, end, is_model)` (unidip.py:98-158): `dat = self.dat[start:end]` and `(_, pval, modidx) = diptst(dat, False, 100)`.
  * `pval is None` returns `[]`.
  * `pval > 0.05`: return `[(start, end-1)]` if `is_model`, else `[(start+w0, start+w1)]` with `w = _get_full_interval((start, end))` (mirror test, unidip.py:160-227).
  * Otherwise recurse:
    * `rm = _unidip(start+modidx[0], start+modidx[1], True)`; `subd = [(t0-start, t1-start+1)]`; `l = min(subd+[modidx], key=second)`; `h = max(subd+[modidx])` (tuple order);
    * if `diptst(dat[:l[1]]).pval < 0.05` (and not None), add `_unidip(start, start+l[0], False)`;
    * if `diptst(dat[h[0]:]).pval < 0.05`, add `_unidip(start+h[1], end, False)`.
* `diptst(dat, is_hist, numt=100)` (dip.py:47-64):
  * `d, (_, idxs, left, _, right, _) = dip_fn(dat)`.
  * `unifs = np.random.uniform(size=numt*len(idxs)).reshape(numt, len(idxs))`. This is the **global, never seeded numpy RNG** (nondeterministic, §5.1).
  * `unif_dips = apply_along_axis(dip_fn, 1, unifs, False, True)`.
  * `pval = None if unif_dips.sum() == 0 else (count(d < unif_dips) + 1)/(numt + 1)`.
  * Returns `(d, pval, (len(left)-1, len(idxs)-len(right)))`.
* `dip_fn` (dip.py:66-135), Hartigan dip after tatome/dip_test:
  * `idxs` = sorted **unique values** (`collections.Counter` keys), `histogram` = their counts.
  * If `len(idxs) <= 4` or all values are equal, return `d=0`, `left=[]`, `right=[1]`.
  * Otherwise `cdf = cumsum(hist, float) / cdf[-1]`, then the GCM/LCM loop (dip.py:93-135) with `_gcm_`/`_lcm_`/`_touch_diffs_` (dip.py:23-45). Transliterate it as-is. Note that `modidx` is in **unique-value index space**, yet `_unidip` uses it as an offset into the full (duplicated) sorted data. Uniform null samples have size `len(idxs)` (the number of unique values), not n.
* Practical behaviour: with large het-site counts, the DP histogram's dip is far below the n≈#unique-values uniform null, so `pval ≈ 1` and the result is `[(0, n-1)]` regardless of the RNG draw.

### 3.2 Coverage loading, masking, genome statistics (analysis.py:104-172)
1. For each genome chromosome with coverage (§2.1), load `positions` (int64 ndarray) and `coverage` (float64 ndarray). On the first such chromosome, set `bin_size` and load the metadata (§2.3.2), then run `pre_calculation_check` (bug B1). Apply the N/blacklist mask (§2.3.2). Store the arrays in `tmp_cov[chrom]` and `tmp_pos[chrom]`. **These are the same array objects that are mutated in place later.**
2. `V = concat(tmp_cov[c] for c in insertion order)` as a python list, turned into a float64 array. It covers **all genome chromosomes with coverage, not just `--only-chr`**, and is taken after masking but before outlier removal.
   * `genome_wide_median = float(np.nanmedian(V))`;
   * `genome_wide_mean = float(np.nanmean(V))`;
   * `genome_std = np.nanstd(V)` (ddof 0).

   The NaN semantics (§6) matter: nanmean/nanstd sum the **full-length** array with NaN replaced by 0.0 (numpy pairwise), then divide by the non-NaN count.
3. `--sample-coverage-overwrite x` sets `genome_wide_median = x` (the mean is not changed). `normalization_value = genome_wide_median` (analysis.py:168-171).
4. If SNV data is present, log a warning when `!(norm - std < snv_med < norm + std)` (no effect).

### 3.3 Per-chromosome normalization (analysis.py:186-192, 221-291)
For `c` in `only_chr_list` order: take `cov = tmp_cov[c]` (KeyError if absent) and `pos = tmp_pos[c]`.
* `norm = snv_med` if SNV was given, else `genome_wide_median`. **`--sample-coverage-overwrite` is ignored when `--snv` is given.**
* Outlier removal (unless `--disable-max-coverage`): for each bin, `if cov > genome_wide_mean + (max_std_outlier_rm * genome_std): cov = NaN` (in place; NaN compares false and stays NaN) (analysis.py:239-243).
* If `nan_count < cov.size` (and `c` is in chr_lengths, always true):
  * `normalized_cov = cov / norm`;
  * `normalized_cov_ploidy = normalized_cov * chr_ploidy[c]` (IEEE elementwise; NaN propagates);
  * store `coverage_analysis[c] = {cov_data, statistics, norm_statistics}`.

  Otherwise the chromosome is silently excluded from `coverage_analysis`. The statistics objects (mean, std, median, min, max) are only used by plots and logs.
* `self.normalization_value` ends up as the last value assigned. It is written to `.spc` (analysis.py:812).

### 3.4 DEL/DUP thresholds (`CNVMetrics`, cnv_metrics.py; analysis.py:200-212)
`CNVMetrics(coverage_analysis, cnv_calls=self.cnv_calls_list (empty dict, same object), exclusion_zones=metadata, hashname=sample_id, ploidy=--ploidy, is_cancer)`.

**Random coverage sample S** (identical every time it is recomputed; see Q-EXCL):
1. `df`: one row per bin of every `coverage_analysis` chromosome, in dict order then position order, with columns `index, chr, position, coverage = normalized_cov_ploidy` (explode, cnv_metrics.py:42-48).
2. Replace ±inf by NaN and drop NaN coverage rows, giving `df_no_nans` with **M** rows (cnv_metrics.py:211-212).
3. The blacklist and CNV exclusion logic (cnv_metrics.py:215-223, 255-275) has **no effect on the selected rows**, because the boolean mask is taken from the unfiltered copy, which is all False (Q-EXCL, §7). It can still crash (B5). The sampled population is all of `df_no_nans`.
4. `random_sample_selection` (cnv_metrics.py:118-135):
   * `random.seed(hashlib.md5(sample_id.encode()).hexdigest())` (the 32-char lowercase hex **string**);
   * `last = M - 1`; `k = int(last * 0.1)`;
   * `idx = sorted(random.sample(range(0, last), k))`. The **last row (M-1) is never selectable**. CPython's exact algorithm is in §6.12.
   * `S = df_no_nans.iloc[idx]["coverage"]` (float64 Series, row order).

**Thresholds** (`calculate_del_dup_borders_quantile(0.1, 0.9)`, cnv_metrics.py:411-435):
```
mean = float(np.mean(S))          # pandas nanmean (pairwise sum / count)
sd   = float(np.std(S))           # pandas nanstd, ddof=0
C    = S[S < mean + 2*sd]         # strict <
DEL  = C.quantile(0.1)            # pandas -> np.percentile(C, 10.0, linear)
DUP  = C.quantile(0.9)
if not cancer:  DEL = ploidy-1 if DEL < ploidy-1 else DEL          # cnv_metrics.py:439-450
                DUP = ploidy+1 if DUP > ploidy+1 else DUP          # (python int replaces float)
```
`del_threshold` and `dup_threshold` are stored. They also become `lower_2n_threshold`/`upper_2n_threshold` (analysis.py:211-212).

**Strict thresholds** (`calculate_del_dup_borders_sd(1.85)`, cnv_metrics.py:387-406), used by clean_by_threshold: `strict_DEL = mean - 1.85*sd`, `strict_DUP = mean + 1.85*sd` (`sd = np.std(S)`, unclipped), clamped the same way.

### 3.5 Raw candidate segmentation (call_cnv_coverage.py:20-89; analysis.py:311-331)
For each `c` in `coverage_analysis` order, scan `zip(normalized_cov_ploidy, positions)` with state `run_start=0`, `run_current_pos=0`, `cur_type=""`, `pos_list=[]`, `cov_list=[]`:
```
for (v, p):
  if isnan(v): continue                                   # NaN bins update nothing
  if v > DUP or v < DEL:
      t = "DUP" if v > DUP else "DEL"
      if run_start == 0:                                  # (A) "start run": p NOT appended
          cur_type = t; run_start = p
      elif p - run_current_pos < bin_size + 1 and cur_type == t:
          pos_list.append(int(p)); cov_list.append(float(v))
      else:
          emit(pos_list, cov_list, cur_type)
          pos_list=[int(p)]; cov_list=[float(v)]; run_start = p; cur_type = t
  else:
      emit(pos_list, cov_list, cur_type)
      pos_list=[]; cov_list=[]; run_start = p; cur_type = ""
  run_current_pos = p
emit(pos_list, cov_list, cur_type)                        # final
emit(...): if len(pos_list) >= 10 -> new CNVCandidate          (call_cnv_coverage.py:13-18)
```
Quirk Q-RUNSTART: branch (A) fires while `run_start == 0`, i.e. until the first non-NaN bin at a position > 0 has been processed. In that state the first out-of-threshold bin is **dropped**. A gap of NaN bins breaks a run because `run_current_pos` is the last *non-NaN* bin.

`CNVCandidate.push_candidates` (cnv_candidate.py:37-50) sets:
* `start = pos[0]`, `end = pos[-1]` (bin **starts**), `size = end - start`, `size_kb = int(size / bin_size)`;
* `median_cov_norm = np.nanmedian(cov)`;
* `cn_status = int(round(median, 0))` (numpy rint, ties to even). The `== np.nan` check is always false. `median == inf` gives `(2, 2.0)` (cnv_candidate.py:99-106);
* `median_raw_cov = np.nanmedian(median*1.0)` (= median);
* `gt` from cn_status (cnv_candidate.py:80-88): 0 or 4 → `1/1`, 1 or 3 → `0/1`, 2 → `0/0`, other → `./.`;
* `type = cur_type`, `sample_origin = sample_id`, `bin_size = int(bin_size)`, `normalization_value = 1.0`.

**cn_status and gt are never recomputed afterwards** for CNVs (grep confirms only cnv_candidate.py:46-48). The final VCF `CN`/`GT` come from the raw candidate object that survives merging.

IDs (analysis.py:324-327): `CNV_ID.n_id_generator(existing_cnv_ids, n=len(cands))` produces, per candidate, 8 × `random.choice("A..Z0..9")`, redrawn if the id is already in `existing` (cnv_id.py:7-36). `cand.id = f"Spectre.{type}.{id}"`. The Python `random` stream continues from the state left by the S sampling of §3.4 (§5.1). `existing_cnv_ids` receives each id twice (cnv_id.py:33 and analysis.py:325), which is harmless.

`raw_cnv_calls_list = cnv_calls_list.copy()` (analysis.py:331) is a **shallow** dict copy: the same list and candidate objects. Later in-place mutation of those objects (merge accumulation, cnv_fit, statistics) is visible in `.spc` `raw_cnvs` (Q-ALIAS).

### 3.6 Merge (analysis.py:361-483), per chromosome
If `len(cands) > 1`:
1. Pass A: `n, M = merge_pass(cands)`; repeat `merge_pass(M)` while `n > 0`.
2. Pass B: `merge_pass(M, allow_over_blacklist=True)` repeated until `n == 0` (always runs at least once).
3. Return `M` if non-empty, else the input (analysis.py:385). With ≤1 candidate the input is returned unchanged.

`merge_pass(c)` (analysis.py:397-483). Pairs are disjoint within a pass:
```
i = 0; out = []; acc = c[0]; last_merged = False; n = 0
while i < len(c) - 1:
    nxt = c[i+1]                                       # note acc is c[i] at loop top
    same_type   = acc.type == nxt.type
    same_status = abs(acc.median_cov_norm - nxt.median_cov_norm) <= cov_diff_threshold   # 0.8
    m   = np.mean([acc.size, nxt.size])                # (0+a+b)/2 in float64
    pct = 40 * np.e ** (-m / 1000000) + 0.1            # pow(e, x) (libm pow, NOT exp)
    thr = m * pct / 100
    dist = max(thr, 10000) + bin_size
    dist_ok = acc.end + dist >= nxt.start
    bl = (chrom in metadata) and any(acc.end <= int(v) <= nxt.start for (s,e) in metadata[chrom] for v in (s,e))
    if allow_over_blacklist: bl = False
    outside = not (DEL <= acc.median_cov_norm <= DUP)  # quantile thresholds
    if dist_ok and same_type and same_status and not bl and outside:
        n += 1
        acc.add_candidates(nxt.pos, nxt.cov, nxt.id, nxt.merged_sample_references)
        out.append(acc); i += 1
        if i >= len(c) - 1: last_merged = True
    else:
        acc.median_coverage_candidates_merged(); out.append(acc)
    i += 1
    if i < len(c): acc = c[i]
if not last_merged: out.append(acc)
```
`add_candidates` (cnv_candidate.py:52-68):
* `pos += other.pos`, `cov += other.cov` (lists; gaps are not filled here);
* recompute start, end, size, size_kb;
* `median_cov_norm = np.nanmedian(cov)`; `median_raw_cov = median_cov_norm * 1.0`;
* `merged_sample_references |= {other.id} | other.refs`.

Type, cn_status and gt are kept from `acc`. The blacklist test checks whether **any coordinate** of any metadata interval (MDR int or blacklist string) lies in `[acc.end, nxt.start]`.

### 3.7 Scaffold (gap filling), analysis.py:666-715
Runs only if the merged list has ≥2 candidates (analysis.py:344-345). For each candidate:
* `d = np.diff(pos)`. For each index `k` where `d[k] != bin_size`, in **descending** order:
  * `gs`/`ge` = the indices in chromosome `positions` equal to `pos[k]` and `pos[k+1]`;
  * `gap = normalized_cov_ploidy[gs:ge]` (includes the candidate's own bin at `pos[k]`; excludes `pos[k+1]`); `m = np.nanmedian(gap)`;
  * if **no NaN** in `gap`: `cov = cov[:k] + [m]*len(gap) + cov[k+1:]` and `pos = pos[:k] + positions[gs:ge] + pos[k+1:]` (the own bin at `k` is **replaced** by `m`).
* Then `median_cov_norm = np.nanmedian(cov)` (0.0 if empty) (cnv_candidate.py:90-94). `median_raw_cov` is **not** updated.

### 3.8 Breakpoint refinement: grow and trim (`CNV_Post_Processing.cnv_fit`, cnv_post_processing.py)
`trim_merged_candidates` (analysis.py:498-511) runs on every non-empty list. `df_cov = DataFrame(position=positions, ploidy_cov=normalized_cov_ploidy)` covers the whole chromosome, NaNs included. For each candidate `X` (mutated **in place**):

**Grow** (cnv_post_processing.py:43-101):
```
n   = max(100*bin_size, X.size*0.25)
cnv = rows with X.start <= pos <= X.end
bef = rows with X.start - n <= pos < X.start
aft = rows with X.end < pos <= X.end + n
E   = concat(bef, cnv, aft)                  # position order
v   = E.ploidy_cov (ndarray)
E.slope_1/5/10 = slope(v, 1/5/10)
slope(v, w): d = zeros(len(v)); for i in range(0, len(v), w):
                 d[i:i+w] = v[i] - np.mean(v[i:i+w+1])     # NaN-propagating; pairwise mean
cnv_sd   = np.std(cnv.ploidy_cov.values) * 3                # ddof 0, NaN-propagating
cnv_mean = np.mean(cnv.ploidy_cov.values)
E.is_in_sd = 1 if cnv_mean - cnv_sd < x < cnv_mean + cnv_sd else 0   # NaN -> 0
up = E.slope_10.mean() + E.slope_10.std()    # pandas, skipna, std ddof=1
lo = E.slope_10.mean() - E.slope_10.std()
ok(r) = r.is_in_sd == 0 and lo < r.slope_10 < up
new_end   = first  pos of rows with ok and pos > X.end    else X.end
new_start = last   pos of rows with ok and pos < X.start  else X.start
Ext = rows of E with new_start <= pos <= new_end
```
If the CNV span contains **any** NaN bin, `cnv_mean` and `cnv_sd` are NaN and every `is_in_sd` is 0 (Q-FITNAN). In that case `new_end` and `new_start` are the nearest bins with a moderate slope.

**Trim** (cnv_post_processing.py:115-142):
* `Ext.slope_diff = abs(Ext.slope_5.diff())`, the difference to the previous **row** of Ext (first row NaN).
* `R = Ext[X.end < pos < new_end]`; `new_end_trim` = **first** pos in R with `slope_diff == R.slope_diff.max()` (skipna), else `new_end`.
* `L = Ext[new_start < pos < X.start]`; `new_start_trim` = **last** such pos with `slope_diff == L.slope_diff.max()`, else `new_start`.

So a candidate can only grow or stay the same.

**Overwrite** (cnv_post_processing.py:144-168): `rows = Ext[new_start_trim <= pos <= new_end_trim]`, then:
* `X.pos = rows.position.tolist()` (contiguous bins, NaN bins included);
* `X.cov = rows.ploidy_cov.tolist()` (raw values; scaffold values are discarded);
* `X.start = pos[0]`, `X.end = pos[-1]`, `X.size = abs(end - start)`.

`median_cov_norm`, `median_raw_cov`, `cn_status`, `gt`, `size_kb` and `type` are **not** updated.

### 3.9 Overlap resolution and final filters
**clean_by_overlap** (analysis.py:513-611) is transliterated exactly, including its index bugs:
```
remove = []; i = 0; off = 0; N = len(L)
while i + off < N:
    cur = L[i+off]; stop = False
    while i + off < N - 1 and not stop:
        nxt = L[i+off+1]
        if not (nxt.start < cur.start < nxt.end or nxt.start < cur.end < nxt.end): break
        cur_full = cur.start <= nxt.start and nxt.end <= cur.end    # unreachable given the test above
        nxt_full = nxt.start <= cur.start and cur.end <= nxt.end
        if cur_full: remove.append(i+off+1); off += 1
        if nxt_full: remove.append(i); stop = True                  # BUG: i, not i+off
        else:
            master, slave, rm = cur, nxt, i+off+1
            if cur.size < nxt.size: master, slave, rm = nxt, cur, i # BUG: i, not i+off
            if master.end >= slave.end:   # "left overlap"
                k = slave.pos.index(master.start)
                master.add_candidates(slave.pos[:k], slave.cov[:k], slave.id, slave.refs, push_front=True)
            else:
                k = slave.pos.index(master.end) + 1
                master.add_candidates(slave.pos[k:], slave.cov[k:], slave.id, slave.refs, push_front=False)
            remove.append(rm); off += 1
    i += 1
return [L[j] for j in range(N) if j not in remove]
```
`list.index` raises ValueError if the position is absent (it is normally present, since all positions lie on the shared bin grid). `add_candidates` recomputes `median_cov_norm` and `median_raw_cov` over the merged `cov` (NaNs ignored by nanmedian).

**clean_by_threshold** (analysis.py:631-656). Thresholds are strict_DEL/strict_DUP (§3.4). Keep candidate `X` iff:
```
mu = np.mean(X.cov)                      # NaN-propagating (cov may contain NaN after cnv_fit)
n_out = count(cov < sDEL) + count(cov > sDUP)
keep = not (sDEL <= mu <= sDUP) and n_out > len(cov) * qc_noise_allow
```
A NaN `mu` passes the first condition.

**clean_by_size** (analysis.py:613-622): keep iff `size >= min_cnv_len * 1.0`.

The resulting list becomes `cnv_calls_list[c]` (analysis.py:352). Its order is the merge order (not re-sorted).

### 3.10 Quality score: z-test (cnv_metrics.py:452-484, 162-199)
`evaluate_cnvs(refined)` first calls `__recalculate_exlustion_zone`, which reseeds and redraws an **identical S** (cnv_metrics.py:460-462, 281-283). For every candidate in `cnv_calls[c]`:
* unrefined (first call): `statistics["z-score"] = {"statistics": None, "score": 0, "pvalue": 0, "sample_score": 0}`;
* refined:
  ```
  x  = array(cov); x = x[~isnan(x)]; x = x[~isinf(x)]
  r  = S.replace(±inf, NaN).dropna()
  mean = np.mean(x); sd = np.std(x)                        # ndarray, ddof 0
  z = (np.mean(r) - mean) / (sd / np.sqrt(len(r)))         # np.float64: /0 -> ±inf or nan, no exception
  zc = min(z, 10) if z > 1 else max(-10, z)                 # py min/max: NaN -> -10 (int); z>10 -> 10 (int)
  p  = scipy.stats.norm.sf(abs(zc)) * 2
  gq = -np.log10(p) * 10
  sample_score = min(60, int(gq))                           # int() truncation
  statistics["z-score"] = {"statistics": None, "score": zc, "pvalue": p, "sample_score": sample_score}
  ```
  Since `len(r) ≈ 0.1·M` (~3e5), `|z|` almost always exceeds 10, giving `p = 2·sf(10) ≈ 1.52e-23` and GQ = 60. GQ < 60 only when `|zc| <~ 4.89`. `norm.sf(x)` = Cephes `ndtr(-x)` (§6.10).

### 3.11 SNV allele-frequency concordance (call_cnv_AF.py:75-126). Computed but **does not affect outputs**
For each refined CNV:
* drop it if `cn_status == 2 and --ploidy == 2`;
* else `af = np.nanmean([AF of records in f"{c}:{start}-{end}" whose AF is not a tuple])` (bin starts used as 1-based coordinates; empty list gives NaN). A `None` AF crashes.
* concordance:
  * ploidy != 2: true;
  * `int(cn) > 6`: true;
  * else `expected = {0:0, 1:0, 3:1/3, 4:[1/4,2/4], 5:[1/5,2/5], 6:[1/6,2/6,3/6]}[cn]`, with `thr = e*0.3` (lists, or scalar with cn > 1) or `0.1` (cn 0, 1). Keep iff any `abs(af - e) <= thr`.

The result goes to `cnv_calls_list_af_filtered` (spectreCNV.py:108-109). That variable is only tested for `is not None` (analysis.py:799); the BED/VCF use the unfiltered `cnv_calls_list` (analysis.py:724). The fetch side effects and the crash paths remain.

### 3.12 SNFJ breakpoint support (breakpoint_analysis.py:27-76), only with `--snfj`
* Load JSON (gzip if the path ends with `.gz`). A missing file returns None and then crashes.
* Matching chromosomes are `set(cnv keys) & set(snfj keys)`. If there are none, warn and skip, but `snfj_breakpoints` is still set True (analysis.py:358), so `SVSUPPORT=FALSE` is written for all records.
* All DEL entries of **all matching chromosomes are pooled** into one DataFrame; the same for DUP (Q-SNFJ). For each refined CNV, with `off = 5*bin_size`: `sv_support = any(row in pool[type] with |svlen| in [size-off, size+off] and pos in [start-off, start+off] and end in [end-off, end+off])`. Non-DEL types use the DUP pool. A matching chrom without a `DEL`/`DUP` key, or with an empty list (no `svlen` column), raises KeyError.

### 3.13 Loss of heterozygosity (snv_analysis.py:255-507; loh_candidate.py)
**Candidate scan** (`candidate_loh`, snv_analysis.py:255-342). For each `(chro, info)` in **VCF header contig order** with `info.length > 1e6` (hard-coded, not `--min-chr-len`) and `chro in only_chr_list`, start with `cand = new LoHCandidate()`, `not_loh = 0`, `loh_per = 0`, `sites = 0`, then iterate `fetch(chro)`:
```
if qual > 5:
    if AF is tuple: continue
    dp, gt, af = int(DP), GT, float(AF)
    if gt == (1,1):
        loh_per += 1; sites += 1
        if not cand.is_init: cand.start_candidate_list(contig, pos, dp, af, float(qual))
        else:
            if 2 > not_loh > 0: (log counter only)
            not_loh = 0
            cand.push_candidates(pos, dp, af, float(qual))
    else:
        not_loh += 1
        if not_loh >= 2:
            if len(cand.pos) >= 10 and loh_per >= sites * 0.70: SAVE(cand, gt_alt_prop = loh_per/sites)
            sites = loh_per = not_loh = 0; cand = new LoHCandidate()
        else:
            sites += 1
# end of contig:
if len(cand.pos) > 10: SAVE(cand, gt_alt_prop = loh_per/sites)      # strict >, no 0.70 test
SAVE: cand.set_copy_number(cn_med, ploidy=2); cand.id = "Spectre.LOH." + n_id_generator(existing_cnv_ids)
```
Records with `qual <= 5` are ignored entirely. Counters are not reset when a candidate starts, so a preceding isolated non-alt site counts in `sites`.

* `start_candidate_list`: `pos=[p]`, `start=end=p` (1-based VCF POS), `size=1`, `size_kb=0.0`, `cov=[dp]`, `aflist=[af]`, `quallist=[q]`, `number_snvs=1` (loh_candidate.py:29-42).
* `push_candidates`: append to the lists; `end=p`; `size=end-start`; `size_kb=int(size/1000)`; `number_snvs += 1` (loh_candidate.py:45-57).
* `set_copy_number(cn, ploidy)` (loh_candidate.py:116-120):
  * `median_cov = np.nanmedian(cov)`;
  * `copy_number = round(median_cov/cn*ploidy, 2)` (numpy round; §6.8);
  * `median_cov_norm = copy_number`;
  * `cn_status = int(round(median_cov/cn*ploidy, 0))`.

**Neighbour merge** (`candidate_loh_neighboring`, snv_analysis.py:344-447). Per contig, in dict order, run up to 500 cycles of `merge_candidates` until a cycle performs no merge:
```
merge_dist(a, b): m = np.mean([a, b]); thr = m * (1/(log10(m)**1.75) + (log10(1e6)/log10(m))/100)
for k, A in enumerate(L):
    if k is last: if combine empty: out.append(A) elif A.id in combine_ids: out.append(COMBINE(combine)); break
    B = L[k+1]
    near = (B.start - A.end) <= merge_dist(A.size, B.size)
    prop = compute_gt_alt_prop(f"{chrom}:{A.start}-{B.end}") if near else 0.0
    if near and prop >= 0.85: combine = combine or [A]; combine.append(B)  (ids likewise); merged = True
    else:
        if combine empty: out.append(A)
        elif len(combine)==1: out.append(combine.pop())           # unreachable
        elif A.id in combine_ids: out.append(COMBINE(combine))
        combine = []
```
* `compute_gt_alt_prop(region)` = `#records with GT==(1,1) and qual>5` / `#all records` in the region (**all** records, including low-QUAL, tuple-AF and indels overlapping the region) (snv_analysis.py:560-568).
* `COMBINE` (loh_candidate.py:59-93): new candidate with `id`, `chrom` and `start` of the first; pos/cov/aflist/quallist concatenated; `number_snvs` summed; `end` of the last; `size`, `size_kb`; `set_copy_number`; `gt_alt_prop` recomputed over `chrom:start-end` as above.

**Filter** (snv_analysis.py:449-507). Applies to the merged list (if non-empty). For each candidate, with `perkb = int(size_kb) * loh_min_snv_perkb`:
* `size > loh_min_region_size`, else `filter = f"LoH region size kb: {size_kb}|{number_snvs}|{perkb}/kb"`;
* `number_snvs > perkb` **or** `number_snvs > loh_min_snv_total`, else `filter = f"Number of SVs: {n}|{perkb}/kb|{size_kb}"`;
* `gt_alt_prop >= 0.85`, else `filter = f"ALT GT propotion: {gt_alt_prop}|0.85"`;
* pass:
  * `loh_pass=True`, `filter="PASS"`, `gt="1/1"`;
  * `het_score = np.round(np.nanmean(aflist), 2)`, `quality = np.round(np.nanmean(quallist), 2)`;
  * `statistics["z-score"]["score"] = ["sample_score"] = min(60, 60/(1+math.exp(-(size/1e6 - 1)))) * np.nanmean(aflist)` (snv_analysis.py:245-252).

`loh_pass_and_non_pass` (snv_analysis.py:528-543) builds dicts keyed by `only_chr_list`, in list order, with candidates in list order. **Pass** goes to `snv_loh` (merged into the outputs); **fail** goes to `snv_loh_raw` (`.spc` only).

### 3.14 CNV + LoH merge for output (loh_candidate.py:152-246; analysis.py:717-730)
`chromosome_names_out = coverage_analysis.keys()`, or `only_chr_list` if that is empty.
* No SNV (`loh is None`): `merged = cnv_calls_list`.
* `len(cnv) == 0` (only with --loh-only): `merged = loh`.
* Otherwise, per chromosome, a two-pointer merge by `start`: when both lists remain, take the CNV if `loh.start >= cnv.start`, else the LoH (ties put the CNV first). CNV starts are 0-based bin starts; LoH starts are 1-based SNV positions. A failsafe counter `nL + nC + 1` cannot trip in practice; if it did, the code falls back to CNV only.

---

## 4. Outputs

All paths are under `out_dir`.

| file | writer | notes |
|---|---|---|
| `{sample}.spc.gz` | analysis.py:784-840, outputWriter.py:245-278 | gzip JSON; written first |
| `{sample}_cnv.bed.gz` + `.tbi` | analysis.py:733-744, outputWriter.py:13-39 | |
| `{sample}.vcf.gz` + `.tbi` | analysis.py:746-751, outputWriter.py:69-242 | **name has no `_cnv`** (`self.output_vcf` from spectreCNV.py:23 is never used) |
| `img/{sample}_plot_cnv_{chrom}.png` | analysis.py:762-775, plot.py:35-86 | one per `coverage_analysis` chrom; skippable |
| `metadata.mdr` (or `--metadata` path) | §2.3.1 | only if missing |
| `tmp/` dir | spectre.py:306-307 | empty (the `snv_af_bed` path is set but never written, spectreCNV.py:105-106) |
| SNV `.vcf.gz`/`.tbi`, coverage `.csi` | §2.4, §1.3 | input-side side effects |
| dev only: `debug/` | `raw-cnvs{sample}.bed.gz` (analysis.py:742-744), `cnv_{sample}_byCoverage_chr{c}.csv` (analysis.py:844-855), `het_dump.txt.gz`, `loh_dump.tsv.gz` (snv_analysis.py:130-135, 509-521), `genome_coverage_distribution-{sample}.png`, `cnv_metrics_random_samples.csv.gz` (written to `{out_dir}/debug/`, cnv_metrics.py:236-240) | |

### 4.1 BED (outputWriter.py:17-39)
For `c` in `chromosome_names_out` order: `sorted(merged[c], key=start)` (python stable sort). Each line is
`"\t".join([chrom, str(start), str(end), type, str(size), str(round(median_cov_norm, 2))]) + "\n"`.
* No header. `start` and `end` are the first/last **bin starts** (CNV) or 1-based SNV positions (LoH).
* `median_cov_norm` is an np.float64, so the value is numpy-rounded (§6.8) and printed with python/numpy shortest repr (§6.9), e.g. `1.0`, `0.53`, `nan`.
  * CNV: the median from the last of push / merge `add_candidates` / scaffold recompute / overlap `add_candidates`. It is **not** the post-`cnv_fit` coverage (Q-STALE).
  * LoH: `copy_number`, rounded a second time.
* The plain file `{sample}_cnv.bed` is written, then `pysam.tabix_compress` (BGZF) and `tabix_index(preset="bed")` produce `.tbi`, and the plain file is removed if the `.gz` is non-empty.

### 4.2 VCF (outputWriter.py:89-242)
**Pre-sort text** written to `{sample}.vcf`. Header lines, in order:
```
##fileformat=VCFv4.2
##FILTER=<ID=PASS,Description="All filters passed">
##source=Spectre
##contig=<ID={c},length={len}>            (one per genome_info["chromosomes"], FASTA order; ALL kept contigs)
##ALT=<ID=DEL,Description="Deletion">
##ALT=<ID=DUP,Description="Duplications">
##ALT=<ID=LOH,Description="Loss of Heterozygosity,">
##FILTER=<ID=UNRESOLVED,Description="An insertion that is longer than the read and thus we cannot predict the full size.">
##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the structural variant">
##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Length of the SV">
##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of copy number variant">
##INFO=<ID=CN,Number=1,Type=Integer,Description="Estimated copy number status">
##INFO=<ID=SVSUPPORT,Number=1,Type=String,Description="Indicator if a SV support was found in a provided SNFJ file">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=HO,Number=2,Type=Float,Description="Homozygosity proportion">
##FORMAT=<ID=GQ,Number=1,Type=Float,Description="Genotype quality">
##FORMAT=<ID=CD,Number=2,Type=Float,Description="Coverage read depth">
##FORMAT=<ID=ID,Number=1,Type=String,Description="Population ID of supporting CNV calls">
##Spectre_sample={sample}
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	{sample}
```
(outputWriter.py:96-137). With zero kept contigs the contig entry is an empty string, which produces an empty header line.

**Records** (outputWriter.py:147-208): `chromosome_names_out` order, then merged list order. Fields are tab-joined:
| col | CNV | LoH |
|---|---|---|
| CHROM | chromosome | chromosome |
| POS | `start` (0-based bin start, written raw) | first SNV POS (1-based) |
| ID | `Spectre.DEL.XXXXXXXX` / `Spectre.DUP.…` | `Spectre.LOH.XXXXXXXX` |
| REF | `N` | `N` |
| ALT | `<DEL>`/`<DUP>` | `<LOH>` |
| QUAL | `.` (`str(".")`) | `str(np.round(mean(QUAL),2))` |
| FILTER | `.` | `PASS` |
| INFO | `END={end};SVLEN={size};SVTYPE={type};CN={cn_status};` (+`SVSUPPORT=TRUE;`/`FALSE;` if `--snfj`) | same, with `SVTYPE=LOH`, `CN` from LoH cn_status |
| FORMAT | `GT:HO:GQ:CD` | `GT:HO:GQ:CD` |
| sample | `{gt}:{round(0.0,2)}=0.0:{int(gq)}:{round(median_raw_cov,2)}` | `1/1:{round(het_score,2)}:{int(sigmoid)}:nan` |

* `END` is the CNV's last **bin start**; `SVLEN = end - start`. The trailing `;` is present in the pre-sort text.
* CNV `HO` is the python float `0.0`, giving `0.0`. LoH `HO`: np.float64 rounded twice, numpy style.
* `GQ` = `int(statistics["z-score"]["sample_score"])`.
* `CD`: CNV `median_raw_cov` (np.float64, from the last `add_candidates` or push; neither scaffold nor fit updates it) with numpy round; LoH `round(np.nan, 2)` gives `nan`.
* The single-sample branch always applies (`support_cnv_calls` is empty, outputWriter.py:165-168).

**Post-processing** (outputWriter.py:229-242): `pysam.bcftools.sort("-o", f"{tmp}.sort.tmp", tmp)`, then `tabix_compress` into `{sample}.vcf.gz`, then `tabix_index(preset="vcf")`, then the temporaries are deleted. The final bytes are therefore **htslib's re-serialization** (to be verified against a real run; these follow from htslib 1.18 behaviour):
* Records are sorted by (header contig index, i.e. FASTA order; POS; then `strcasecmp` of REF, ALTs), qsort (§5.6).
* bcftools appends `##bcftools_sortVersion=<ver>+htslib-<ver>` and `##bcftools_sortCommand=sort -o <abs tmp>.sort.tmp <abs tmp>; Date=<ctime>` before `#CHROM` (non-reproducible).
* INFO: the trailing `;` is dropped by htslib's INFO parser. Integers are re-printed.
* Float FORMAT/QUAL are parsed to **float32** and printed via htslib `kputd` (≈`%g`, 6 significant digits): `0.0`→`0`, `60`→`60`, `1.23`→`1.23`, `nan`→`nan`. QUAL `.` stays missing.
* `POS 0` (a CNV starting at bin 0) round-trips as `0`.

### 4.3 `.spc.gz` intermediate (analysis.py:784-840; outputWriter.py:257-278)
`json.dump(obj, gzip.open(path,'wt'), indent="\t")` (python json: separators `,` and `: `, `NaN`/`Infinity` allowed, floats via `float.__repr__`, np.float64 included; a TypeError is caught and logged, which can leave a truncated file). The gzip header holds mtime and the filename, so the bytes are non-reproducible. Object:
```
{"metadata": {"source":"spectre","spectre_version":"0.2.alpha","bin_size":int,"min_cn_len":min_cnv_len,
              "normalization_value":float(normalization_value),"sample_id":sample},
 "genome_info": {"chromosomes":[...],"chr_lengths":{...},"chr_ploidy":{...}},
 "raw_cnvs":        {chrom:[candidate.__dict__...]}   # raw dict: aliased, mutated objects (Q-ALIAS)
 "refined_cnvs":    {chrom:[...]}                     # cnv_calls_list (NOT AF-filtered)
 "refined_loh_cnvs":{chrom:[...]}                     # PASS LoH   (both {} if no --snv)
 "raw_loh_cnvs":    {chrom:[...]}                     # FAIL LoH
 "analysis_metrics": {"min_chr_length":100000,"max_std_outlier_rm":int,"mosdepth_cov_genome_chr_diff":float,
     "lower_2n_threshold":DEL,"upper_2n_threshold":DUP,"strict_del_threshold":sDEL,"strict_dup_threshold":sDUP,
     "cov_diff_threshold":float,"dist_proportion":float,"dist_min_overwrite":10000,"candidate_final_threshold":min_cnv_len,
     "genome_mean":float(np.mean(S)),"genome_sd":float(np.std(S)),"genome_var":float(np.var(S))}}
```
`strict_*` fall back to DEL/DUP if `clean_by_threshold` never ran (analysis.py:805-807).

Candidate dicts are `dict(candidate.__dict__)` with sets converted to lists (`merged_sample_references`, order = python set iteration, **hash-randomized**). Key order is the `__init__` order (cnv_candidate.py:6-35): `chromosome, start, end, size, bin_size, pos, cov, normalization_value, type, cn_status, gt, median_cov_norm, median_raw_cov, size_kb, is_init, filter, quality, id, het_score, statistics, support_cnv_calls, sample_origin, merged_sample_references, sv_support, merged_sample`. LoH appends `aflist, quallist, median_cov, copy_number, pretty_size, number_snvs, gt_alt_prop, gt_alt_count, gt_all_count, loh_pass, af` (loh_candidate.py:15-27), and its `statistics` defaults to the string `"nan"` values (loh_candidate.py:13).

### 4.4 Plots
`CNVPlot` (plot.py:35-86): 8×4 in, dpi 300. It shows `normalized_cov_ploidy` vs position (ylim 0–6), CNV bars (DUP `#d73027`, DEL `#1a9850`) on a lower panel, the mean line and the DEL/DUP lines. Plots can be skipped in a port. `cnv_plot` multiplies `norm_statistics.median` in place (analysis.py:768), which has no other effect.

---

## 5. Nondeterminism and platform dependence

1. **Python `random` (MT19937).** It is seeded with `md5(sample_id).hexdigest()` before every S draw (cnv_metrics.py:126-127), so S and all IDs are deterministic. The stream order is:
   1. seed + sample(S);
   2. CNV ids (8 `choice` per raw candidate, chromosome order);
   3. reseed + sample (unrefined eval);
   4. reseed + sample (refined eval);
   5. LoH ids.

   Because step 5 replays the same stream as step 2, the first K LoH draws collide with the K CNV ids and are redrawn. **LoH ids are stream ids K+1, K+2, …** With `--loh-only` no seed is ever set, so ids come from OS entropy (and the run crashes anyway, B2).
2. **numpy global RNG** (`np.random.uniform` in dip.py:53) is never seeded, so the SNV CN-neutral estimate (§3.1) is nondeterministic in borderline dip cases. A port should use a seedable RNG and accept parity only when `pval` is far from 0.05.
3. **Hash randomization** (PYTHONHASHSEED) affects the order of `merged_sample_references` in `.spc`. The SNFJ chromosome intersection order (breakpoint_analysis.py:34) does not change results.
4. **pandas/numpy reduction order.** numpy pairwise summation (§6.1) is used everywhere a float mean or std is taken. `np.nanmean`/`np.nanstd` sum the **full array with NaN→0**, so the result depends on NaN positions (do not compact first). pandas `Series.mean/std/var/quantile` on float64 use numpy pairwise sums **unless bottleneck is installed**, in which case they switch to sequential sums. The S coverage column is `object` dtype after `explode` and becomes float64 only because pandas 2.x `replace` downcasts object blocks (cnv_metrics.py:211, 259). With pandas 3.x it would stay `object`: sums become sequential python-float additions and quantile runs on objects.
5. **libm.** `pow` (merge distance), `log10` (GQ, LoH merge distance), `exp` (LoH sigmoid) and `log` (CPython `random.sample` setsize) come from glibc. A Rust port on Linux/glibc using `f64::powf/log10/exp/ln` matches; other platforms may differ by an ulp. Note `np.e ** x` is `pow`, not `exp`. `scipy.stats.norm.sf` uses Cephes `ndtr`/`erfc` (not libm).
6. **bcftools sort** uses `qsort` on equal keys (same contig, POS, REF, ALT). glibc's qsort was stable (mergesort) before 2.37; later versions may not be. Ties are rare (duplicate starts after overlap cleanup).
7. **Compression bytes.** BGZF output depends on the zlib/libdeflate build of pysam's htslib, and `.spc.gz` embeds mtime. Compare decompressed content.
8. **Multiprocessing** only in `--population` (`Pool.map` preserves order).
9. Dict and list orders are insertion-ordered (python ≥3.7), so they are deterministic.
10. **Matplotlib-only paths** (plots, the dev histogram) can be skipped. They have no feedback into results, except that `cnv_plot` mutates `norm_statistics.median` after all outputs are written.

---

## 6. numpy / pandas / scipy / CPython semantics to reproduce

| # | call | where | exact semantics |
|---|---|---|---|
| 6.1 | `np.sum` / `np.mean` on contiguous float64 | everywhere | **pairwise**: n<8: `s=0.0; s+=a[i]`. 8≤n≤128: `r[0..8]=a[0..8]`; for blocks of 8 `r[j]+=a[i+j]`; `s=((r0+r1)+(r2+r3))+((r4+r5)+(r6+r7))`; then add the remaining `n%8` sequentially. n>128: `n2=n/2; n2-=n2%8; sum(a[:n2])+sum(a[n2:])`. Mean = sum/n. Integer input → float64 cast in buffered chunks of 8192 (exact for int DP sums). |
| 6.2 | `np.nanmean(x)` | analysis.py:167,251,275; snv_analysis.py:104,114,122,127,252,476-477; outputWriter (pop) | NaN→0 over the **full** array, pairwise sum, ÷ non-NaN count. Ints: plain mean. |
| 6.3 | `np.nanstd(x)` (ddof 0) | analysis.py:172,251,275; snv_analysis.py:115,123,128 | `avg=sum(x0)/cnt`; `d=x-avg` with NaN positions set to 0; `sum(d*d)/cnt`; sqrt. |
| 6.4 | `np.std(nd)` / `np.var(nd)` (ddof 0) | cnv_post_processing.py:72; cnv_metrics.py:181; analysis.py:833-835 (Series) | `m=sum(x)/n; sum((x-m)*(x-m))/n`; NaN-propagating for ndarrays. |
| 6.5 | `np.nanmedian` / `np.median` | cnv_candidate.py:47,64,92,101; analysis.py:166,251,276,697; loh_candidate.py:117 | remove NaN; odd: middle element; even: `(lo+hi)/2` (computed as `0+lo+hi`, then `/2`). All-NaN → NaN (warning). |
| 6.6 | `np.percentile(a, q)` linear / `Series.quantile(q)` | snv_analysis.py:139 (q=98); cnv_metrics.py:429-430 (q=0.1, 0.9 → `np.percentile(v, q*100.0)` = 10.0 and 90.0 exactly) | sort; `v=(n-1)*(q/100)`; `i=floor(v)`; `g=v-i`; `a=s[i]`, `b=s[min(i+1,n-1)]` (if `v >= n-1` both are `s[n-1]`); `lerp = b-(b-a)*(1-g)` if `g>=0.5` else `a+(b-a)*g`. |
| 6.7 | `Series.mean()`, `Series.std()` (**ddof=1**) | cnv_post_processing.py:79-80 | skipna: NaN→0 copy, pairwise sum ÷ count; var = `sum(((avg - x)²) with NaN→0)/(count-1)`; sqrt. `np.mean(Series)`/`np.std(Series)`/`np.var(Series)` dispatch to pandas with ddof=0 (cnv_metrics.py:316-323, 393-395, 420-422, 184; analysis.py:831-835). |
| 6.8 | `round(np.float64, d)` / `np.round(x, d)` | cnv_candidate.py:106; loh_candidate.py:118,120; snv_analysis.py:476-477; outputWriter.py:29,166,168 | **numpy**: `rint(x*10^d)/10^d` (ties to even on the scaled product). Differs from python's correctly rounded `round(float, d)`, e.g. 0.075 → numpy 0.08, python 0.07. Python `round(float, 2)` applies only to CNV `het_score` (0.0) and the LoH `median_raw_cov` (nan). `int(round(x, 0))` = rint, then truncation. |
| 6.9 | `str(np.float64)` / `repr(float)` / f-string | BED col 6, VCF QUAL/HO/CD, JSON | shortest round-trip digits; positional iff `1e-4 <= |x| < 1e16` (or 0), keeping at least one fractional digit (`2.0`), else scientific `d[.ddd]e±XX` (≥2 exponent digits). `nan`, `inf`, `-0.0`. |
| 6.10 | `scipy.stats.norm.sf(x)` | cnv_metrics.py:188 | `ndtr(-x)`, Cephes: `t=x/√2`; `|t|<1/√2` → `0.5+0.5·erf(t)`; else `y=0.5·erfc(|t|)`, `1-y` if `t>0`. Port the Cephes `erf`/`erfc` polynomials so `int(-10·log10(p))` matches. |
| 6.11 | `np.diff`, `np.where`, `np.count_nonzero`, `np.isnan`, `np.cumsum` (sequential), `Series.diff()`, `abs`, `Series.max()` (skipna; empty → NaN), `Series.between` (inclusive both), `DataFrame.explode`, `.iloc`/`.loc` | §3.4–3.9, dip.py | as stated in the sections. Note the iloc-with-labels bug (B5). |
| 6.12 | CPython `random.seed(str)` | cnv_metrics.py:127 | `a = int.from_bytes(s.encode() + sha512(s.encode()).digest(), 'big')`; `init_by_array` with the 32-bit little-endian words of `a` (24 words here). |
| 6.13 | `random.sample(range(n), k)` | cnv_metrics.py:133 | `setsize = 21 + (4**ceil(log(3k)/log(4)) if k>5 else 0)` (math.log → libm log). If `n <= setsize`: pool algorithm (`j=randbelow(n-i)`; `res[i]=pool[j]`; `pool[j]=pool[n-i-1]`). Else: set algorithm (`j=randbelow(n)`, redrawn while already selected). `randbelow(n)`: `kb=n.bit_length()`; `r=getrandbits(kb)` until `r<n`; `getrandbits(kb≤32) = genrand_uint32() >> (32-kb)`. For a ~2.9M-bin genome the set algorithm applies. Both branches are reachable (pool when `0.83·4^c < n ≤ 4^c+21`). |
| 6.14 | `random.choice(36 chars)` | cnv_id.py:14 | `randbelow(36)`: `getrandbits(6)` (`genrand>>26`), rejected if ≥36. Alphabet `ABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789`. |
| 6.15 | `math.exp`, `np.log10`, `np.e**x`, `np.sqrt` | snv_analysis.py:252, 351; cnv_metrics.py:184,190; analysis.py:423 | glibc libm (§5.5). |
| 6.16 | `np.random.uniform` + `np.apply_along_axis` | dip.py:53-55 | legacy global `RandomState` (MT19937, `random_sample`: 53-bit doubles `(a>>5, b>>6)`); unseeded. |
| 6.17 | python `min`/`max` with NaN | cnv_metrics.py:186 | `max(-10, nan)` → -10; `min(nan, 10)` is unreachable (the guard `z>1` is false for NaN). |
| 6.18 | `htslib kputd` (float32 → text) | final VCF | special-cases 0 and ±; outside `[1e-4, 999999]` uses `%g`; otherwise a 6-significant-digit fixed formatter with trailing zeros stripped (kstring.c). |

---

## 7. Bugs, quirks, dead code

**Bugs and quirks that affect output (reproduce first, fix behind a flag later):**
* **B1** `pre_calculation_check` returns `(min_cnv_len <= 1e5) AND (min_cnv_len > 1e4)`. Any `--min-cnv-len > 100000` calls `sys.exit(1)` with no error message (analysis.py:86-101, 146-147).
* **B2** `--loh-only` always crashes: `cnv_metrics` is None at analysis.py:805. IDs would also be unseeded.
* **B3** VCF CN, GT and CD are stale (§3.5, §3.8). `cn_status`/`gt` come from the raw run median of the surviving accumulator candidate. `CD` is the median at the last `add_candidates`. BED col 6 is the pre-fit (scaffolded) median. None reflect the final `cov`.
* **B4** The AF-concordance filter result is never used for output (analysis.py:724 vs spectreCNV.py:109).
* **B5** Blacklist exclusion uses `.iloc[labels]` (labels used as positions, cnv_metrics.py:219-220). It flags the wrong rows and raises IndexError when a label ≥ the row count, i.e. a blacklist start `s` with `s-1` on the bin grid whose bin is non-NaN late in the genome.
* **Q-EXCL** Exclusion masks are taken from the unfiltered frame (all False) and reindexed (cnv_metrics.py:223, 275), so blacklist and CNV regions are **never** excluded from S. Every S draw is identical.
* **Q-RUNSTART** The first out-of-threshold bin is dropped while `run_start == 0` (call_cnv_coverage.py:56-58).
* **Q-UNIDIP** Intervals from sorted, unique-value space are applied to the unsorted DP list with an exclusive end (snv_analysis.py:111-123). The normal result is the median of all but the last het DP. The SNV ploidy is forced to 2 (snv_analysis.py:200). `--sample-coverage-overwrite` is ignored when `--snv` is given.
* **Q-MDR** An N run of length `n_size-1` qualifies; a leading run starts at 1 instead of 0; a trailing run is dropped; header split on space only; `;` lines are treated as headers; GC% uses the last contig's length.
* **Q-BL** Blacklist BED starts are not converted from 0-based; the masking range is inclusive at both ends; the last character of each line is chopped; values stay strings (merge dedup never matches MDR ints).
* **Q-COORD** CNV POS/END are 0-based bin starts and END is the *start* of the last bin. LoH POS/END are 1-based. SVLEN = end − start.
* **Q-SNFJ** Breakpoints from all matching chromosomes are pooled.
* **Q-LOH** The in-loop save requires `len(pos) >= 10` plus the 0.70 test; the end-of-contig save requires `> 10` and skips the 0.70 test. The `or` in the SNV-count filter (snv_analysis.py:469) is probably meant to be `and`. `gt_alt_prop` of combined candidates counts every record (any QUAL, indels) in the denominator. The LoH contig filter is a hard-coded 1 Mb on the VCF header length. Missing header contig lengths crash.
* **Q-VAF** DeepVariant `VAF` (Number=A) makes every record a tuple, which ends in a crash.
* **Q-FITNAN** Any NaN inside a CNV disables the `is_in_sd` criterion in grow (§3.8).
* **Q-OVERLAP** `clean_by_overlap` uses `idx` instead of `idx+offset` for removals, and its full-containment branch is unreachable (analysis.py:542-551).
* **Q-ALIAS** `raw_cnvs` shares objects with refined CNVs: the merge accumulators carry merged pos/cov, fitted coordinates and refined statistics.
* **Q-PLOIDY** `--ploidy-chr c:p` sets the multiplier `int(ploidy²/p)`, not the ploidy.
* **Q-SCAF** Scaffolding replaces the candidate's own bin at the gap start with the gap median, and that median includes the bin itself.
* README defaults for the LoH options disagree with the code (§1.2).
* Version strings disagree (§0).

**Dead code in single-sample calling:**
* `MosdepthReader`, `MosdepthSummary` (mosdepthReader.py), `OSUtil`, `NormaldataAnalyser` (dataAnalyzer.py);
* `VCFParser.vcf_pysam/vcf_to_dataframe/get_mosdepth_chromosome_borders`, `VCFtoCandidate` (vcf_parser.py:47-249);
* `CNVMetrics.calculate_del_dup_borders` (the iterative version, cnv_metrics.py:296-385), `get_ks_test` (137-160), `get_lower/upper_border`;
* `clean_by_coverage` (analysis.py:624-629), `coverage_plot` (analysis.py:754-760, which would crash), `Spectre.meta_data_extraction` (spectre.py:228-262), `FastaRef.extract_n_regions_from_report`, `IntermediateFile.convert_genome_info_to_dictionary`;
* `LoHCandidate.set_gt_alt_prop2`, `vcf_line`; `call_cnv_AF.get_CNV_type_from_af`, `check_af_duplication` (which would crash: `np.nan()`);
* `VCFOutput.id_generator`; the population branch of `vcf_result` (outputWriter.py:169-206); `Breakpoints.__row_constrains`;
* the `dev` subcommand; `--threads` without `--population`; CLI `-03/-04/-06`; `UniDip.plot`, `test_unidip`.

---

## 8. Population mode (summary only)

* **Entry points.** `CNVCaller --population` runs each sample with `Pool(threads).map` (per-sample coverage path made absolute, spectre.py:331-346). If more than one `.spc` results, `population_exe("population_file", paths, out_dir, reference, as_dev)` is called; note `as_dev` lands in the positional `disable_quality_filter` slot (spectre.py:349-350 vs 368-369). `Population --candidates a.spc b.spc ...` calls `population_exe()`. `--reciprocal-overlap` is parsed but ignored: 0.8 is hard-coded (spectre.py:397). The path for non-`.spc` input calls `make_genome_info(reference)` with a bad signature and would crash (spectre.py:380-382).
* **Loading** (spectreCNVPopulation.py:51-167). JSON `.spc` is converted back into `CNVCandidate`/`LoHCandidate` objects by setting `__dict__` keys, with `sample_origin` set to the sample id (deduplicated with a `_k` suffix). `refined_cnvs + refined_loh_cnvs` become `final_candidates`; with `--disable-quality-filter`, `raw_cnvs + raw_loh_cnvs` become `raw_candidates`. Bin sizes must match.
* **Similarity merge** (spectreCNVPopulation.py:356-471). Repeated while it returns more than 1 merge group: per (type, chromosome, cn_status), candidates whose start and end are both within `20·bin_size` of each other are collapsed onto the smallest dataframe index. The superseded ones are marked `merged_sample`, moved into the base candidate's `support_cnv_calls`, and removed; the base `quality` is the max.
* **Support search** (spectreCNVPopulation.py:253-300). For every pair of candidates (any sample, different id) with the same `cn_status` and `type` and reciprocal overlap ≥ 0.8 (`overlap/size` for both), the first is added to `cnv_call_list[chrom]` (a **set**) with support recorded per sample. Merged groups are added back afterwards (spectreCNVPopulation.py:486-508). Optionally, raw candidates are searched as well (spectreCNVPopulation.py:302-335).
* **Output** is `{out}/population_mode_{sample_id}.vcf.gz` via the same `VCFOutput`, in population form:
  * header `##Spectre_population_samples=` plus the `SUPPVEC`/`POPMERGENR` INFO lines;
  * per-sample FORMAT `GT:HO:GQ:CD:ID` with the most common GT, the mean GQ (`np.float16` of the scores), the mean CD and the comma-separated supporting ids, or `./.:0.0:0:0.0:NULL`;
  * INFO appends `;SUPPVEC=0101…;POPMERGENR=n` (outputWriter.py:169-206). This produces a `;;` after the base INFO's trailing `;`.

  Record order before `bcftools sort` follows set iteration (hash-dependent).
