---
name: scanvi-annotation
description: Run OmicsAgent scANVI reference label transfer for scRNA-seq or single-cell spatial data, inspect posterior weights, and audit candidate labels.
---

Use `annotation_method: "scanvi"` with `reference`, `reference_label_key`, and explicit raw count layers. Install the `scanvi` extra. Reject normalized/fractional input counts, duplicate gene IDs, missing labels and inadequate shared genes. Query and reference must share species and compatible gene identifiers; review tissue/label coverage before interpreting results.

Run `omicsagent run --config CONFIG.json`, or the bounded `annotate` tool after marker analysis. The adapter intersects genes, fits scVI then scANVI on a reference/query joint object with origin as batch, and outputs posterior weights, argmax candidates, entropy/margin and saved model. Origin adjustment is not a substitute for donor-aware validation. CPU training defaults to 50 + 30 epochs. The two-epoch demo only checks execution.

Inspect `annotation_weights.csv.gz`, `annotation_gate.csv`, `cells.csv`, and `scanvi_model/`. Softmax probabilities are not calibrated correctness probabilities; scANVI cannot invent a missing reference cell type. Preserve unresolved/unsupported calls. Never treat agreement with TACCO as independent biological proof. SPLIT is opt-in for imaging only, preserves counts and does not rerun clustering on corrected values.

Source adaptation: supplied xenium-scanvi-annotation skill. API: https://docs.scvi-tools.org/en/1.3.0/api/reference/scvi.model.SCANVI.html . Do not repeat private benchmark accuracy claims.
