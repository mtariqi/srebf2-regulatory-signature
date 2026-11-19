# SREBF2 Phase 13–14 reproducibility audit — working checklist

Status: Git feature branch published at `0b11859`; untracked files intentionally excluded.

## Verified in user-provided terminal logs

- [x] Phase 14A.5: 1,102 chromosome-exact pairs, 23 chromosomes in both comparisons.
- [x] Phase 14A.5 SHA-256: `2217e60ed19d244f2f639b5f062c7a26f9d11bf5755a34c645a64cfbe2329270`.
- [x] Phase 14A.6: 2 tests, 3 TSV lines.
- [x] Phase 14A.6 SHA-256: `156bda973a5a03edf4998ef5c00555ba4fedbf2917e2893f882fa54a36ec3778`.
- [x] Git commit `0b11859` pushed to `feature/phase7-cut-and-tag`.
- [x] Git whitespace check passed for staged commit.

## Outstanding before release

- [ ] Capture new SHA-256 of normalized Phase 14A.1 TSV (its pre-normalization checksum is obsolete).
- [ ] Verify Roadmap E117 BED.gz local SHA-256 `f1b30d0fa472c61642b35e2819d130aec0b0359d9bf6abdbdf66418f18271508` again before archival.
- [ ] Add source/reference provenance Markdown, source retrieval commands and licensing information; review reference assembly and liftover assumptions.
- [ ] Document executable dependency versions, environment lockfile or container digest.
- [ ] Version Phase 14A.3 LOCO and Phase 14A.7/14A.8 read-only scripts plus derived summary tables; do not overwrite frozen results.
- [ ] Include targeted unit/integrity checks for matched-pair uniqueness, chromosome agreement, zero-based TSS mapping and counts.
- [ ] Confirm that all state-selection and post-hoc multiplicity caveats appear in the final report.
- [ ] Independently rerun end-to-end workflow on a clean checkout with verified external inputs.

## Reference

`https://egg2.wustl.edu/roadmap/data/byFileType/chromhmmSegmentations/ChmmModels/coreMarks/jointModel/final/E117_15_coreMarks_hg38lift_mnemonics.bed.gz`

Historical HeLa-S3 E117 is not a matched SREBF2 perturbation assay.
