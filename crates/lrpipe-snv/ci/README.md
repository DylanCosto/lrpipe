# Checks

lrpipe's `tools/check.sh` is the check to pass before pushing, run on this machine; GitHub runs it
too, where the GPU and parity checks skip for want of a GPU and model. For lrpipe-snv it runs:

1. `cargo clippy -D warnings` with and without the `gpu` feature (`--no-default-features`).
2. `cargo test --release`, and `snv-infer`'s tests without the `gpu` feature. The GPU tests run
   when the CUDA driver, a GPU and NVRTC are found: NVRTC (`libnvrtc.so.12`) is linked from
   `$NVRTC_ENV/lib` into a directory of its own on `LD_LIBRARY_PATH`, so the rest of that env does
   not shadow other libraries. Where it finds NVRTC and a GPU it sets `SNV_REQUIRE_GPU`, so the
   GPU tests fail rather than skip; kernels that do not compile always fail them.
3. `ci/parity.sh`: `lrpipe snv` on the ONT fixture below, from the BAM and from the CRAM,
   compared with the VCFs Clair3 2.0.3 wrote for the same inputs (`merge_output`, `pileup` and
   `full_alignment`: 467, 2,104 and 603 records), ignoring only the header lines that name the
   command and the reference. It takes the model from `$SNV_MODEL` (this machine's Clair3 env's
   `r1041_e82_400bps_sup_v500` by default), and is skipped, saying so, where there is none. A
   third run, from the BAM with `--gvcf` (as `lrpipe run --gvcf` runs it), must write
   `merge_output.gvcf.gz` with its `.tbi`, holding the records of Clair3 2.0.3's gVCF of the
   same reads (`expected/merge_output.gvcf.gz`: position, alleles, filter, INFO and genotype,
   and a block's GQ, MIN_DP and PL) in order, but for the Clair3 bugs lrpipe-snv fixes
   (`expected/merge_output.gvcf.fixes`, as `diff` shows them; see the README's "Clair3 bugs
   fixed on purpose"): the SNV inside a deletion that Clair3 leaves out, and the blocks of
   chunks without reads, which end before the next chunk's first position and at the contig's
   end.

The CRAM must give exactly the BAM's VCFs. lrpipe-snv's networks sum in another order than
Clair3's libtorch, so a VCF need not be byte-identical to Clair3's: one that is is reported as
"same as the expected output", and otherwise its PASS calls (position, alleles, genotype;
`tools/concordance.py`) must agree with Clair3's on at least 99.9%, with a Jaccard index of at
least 99.9%. On this machine all three VCFs agree 100% but are not byte-identical.

## Fixture

`ci/fixture` holds the 152 reads of HG003 (ONT R10.4.1 e8.2 400 bps, dorado sup v5.0.0) that
overlap chr20:30,200,000-30,240,000, as a BAM and a CRAM, and a chr20 that is GRCh38 from
30,140,001 to 30,310,000 and `N` elsewhere (bgzipped; the contig keeps its length, so coordinates
are GRCh38's). The BAM's `@PG` lines are removed. `tools/make_ci_fixture.sh` rebuilds it and the
expected VCFs. `window.bed` is a 20 kb window of it, and `expected/window/merge_output.vcf`
Clair3 2.0.3's calls with `--bed_fn window.bed` (run by hand), which `lrpipe run --regions`'s
test (`lrpipe/tests/run.rs`) holds lrpipe-snv's PASS calls to.

Check the terms of the reads' source before making the repository public.

The parity check alone, from lrpipe's directory:

```
lrpipe-snv/ci/parity.sh "target/release/lrpipe snv" PATH/TO/r1041_e82_400bps_sup_v500
lrpipe-snv/ci/parity.sh target/release/lrpipe-snv PATH/TO/r1041_e82_400bps_sup_v500
```
