# v0.3.0 validation record

Date: 2026-10-24. Local platform: macOS arm64, Python 3.12.14, Scanpy 1.11.5, AnnData 0.12.19, harmonypy 0.0.10. Optional installed methods: scvi-tools 1.3.3 and TACCO 0.5.0. Exact versions are locked in uv.lock and recorded in run manifests.

## Hierarchical method

The independently generated synthetic reference and query each contain 240 observations and 116 genes, with two L1 lineages and four subtypes. Sixty query observations have deliberately wrong initial L1 labels.

| Run | Rounds | Stop reason | Initially wrong L1 corrected | Final resolved L2 |
|---|---:|---|---:|---:|
| Imaging spatial | 2 | exact_stability | 60 / 60 | 240 / 240 |
| scRNA with Harmony | 2 | exact_stability | 60 / 60 | 240 / 240 |
| Visium domain code path | 2 | exact_stability | 60 / 60 | 240 / 240 |

These are engineered integration fixtures, not biological performance estimates. The Visium fixture is not a mixture/spot benchmark. Real ISD atlas validation is pending; no real ISD data or private case studies are distributed.

## Automated and manual checks

Tests cover old workflows plus marker/atlas conflicts, atlas-only rejection, exact movement threshold, cycle/cap distinctions, stale subtype invalidation, missing markers, constant pools, zero-margin ties, empty reports, cached hierarchy replay without refitting, artifact tampering, raw-count preservation, sample-isolated spatial edges, real Harmony execution, mock-agent hierarchy orchestration and publication scanning. The optional TACCO integration test executes when TACCO is installed. See CI for platform-specific current results; a skipped optional test is not a validation claim.

Ruff and mypy are run on the release source. All project skills are schema-validated. Final-round marker dotplots were inspected visually. The installed package is tested in ordinary non-editable mode as well as the source tree. Editable install paths can be skipped by Python when a local filesystem marks the path file hidden; the tutorial uses --no-editable for portability.

## Limits

No real Claude request was made. Mock-agent tests exercise actual numerical tools. No native R SPLIT parity, real ISD annotation, donor holdout, clinical validation or superiority benchmark is claimed. Earlier optional scANVI two-epoch runs check training and serialization mechanics only. Published history and files are screened for CJK text, common credential formats, private-key blocks, private paths and prohibited local artifacts; credentials are never intentionally included.
