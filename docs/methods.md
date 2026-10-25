# Scientific and engineering contracts

## Expression

Raw counts must be nonnegative, finite and integer-like. Unique observation and feature names are required. QC counts are computed before gene filtering; reported low-gene/high-mitochondrial counts can overlap. Cells must survive the combined mask, and zero libraries after gene filtering cause an explicit error.

For scRNA and Visium, total normalization to 10,000 then log1p is the default. Imaging supports area normalization (counts/area times median area) or log1p without library scaling; total normalization is available with a targeted-panel warning. Mitochondrial prefixes default to human `MT-` and mouse `mt-`; stable Ensembl identifiers need preprocessing or an explicit prefix and review.

PCA uses centered sparse log expression without unit-variance gene scaling; this is an explicit v0.1 choice that avoids densifying the full matrix. For scRNA/Visium, HVGs are selected on a working copy. Imaging uses all retained panel genes. Leiden is igraph, resolution configurable, two iterations. UMAP/graph seeds are recorded. These do not guarantee identical outputs across library versions.

Marker ranking uses full-gene log expression, never scaled PCA input. Singleton clusters are excluded from ranking; insufficient groups produce an empty marker table and warning. Post-clustering p-values are exploratory, not confirmatory biological differential-expression tests. Default broad marker detection is a conservative heuristic, not a trained annotator.

References: [Scanpy rank_genes_groups](https://scanpy.readthedocs.io/en/1.11.x/generated/scanpy.tl.rank_genes_groups.html), [Scanpy Leiden](https://scanpy.readthedocs.io/en/1.11.x/generated/scanpy.tl.leiden.html).

## Spatial

Coordinates must be finite 2D values and units must be explicitly declared. No cross-sample edges are permitted; disconnected/isolated observations are counted. Radius is in the declared units. Independent sections should receive distinct graph sample IDs; this is not automatically the biological donor identifier. Within-sample results are not treated as independent donor-level replication.

Neighborhood enrichment uses symmetric binary radius graphs and sample-specific label permutations. Only the upper triangle of label pairs is emitted to avoid duplicated hypotheses. Same-label `directed_count` counts both directions; between-label values count the row-to-column direction. Two-sided empirical p uses `(1 + exceedances)/(B+1)` around the null mean. Degenerate nulls are undefined. BH adjustment is across all valid hypotheses in each output table. With B=99 the minimum p is 0.01; use more permutations for exploratory resolution, not as a substitute for biological replication.

Moran I uses the binary symmetric graph, permutation of centered gene values within a sample, and a one-sided alternative for positive autocorrelation. Only the configured number of highest-variance genes is tested; BH applies to those emitted tests. This is not a genome-wide SVG screen. No histology segmentation, geometry barriers or covariate-adjusted null is modeled.

Center niches retain fixed center positions and shuffle labels of all neighboring candidates within each sample. Overlapping neighborhoods count incidences. The test asks about labels near those fixed centers, not whether expression-derived selection of centers itself was independent. For formal comparisons use donor-level replicated modeling beyond this version.

Reference API comparison: [Squidpy spatial graph conventions](https://squidpy.readthedocs.io/en/stable/api/squidpy.gr.spatial_neighbors.html). OmicsAgent implements its smaller radius/permutation engine with SciPy; it does not claim byte-for-byte Squidpy equivalence.

## Agent and evidence

Numerical tools are state constrained and usable without a model. The agent's inputs are summaries and fixed schemas. The only external endpoints are Anthropic (when agent mode is requested), Ensembl and Europe PMC (when online evidence is enabled). Retrieval errors are bounded; cache is run-local. Cached evidence retains retrieval timestamps and has no automatic freshness refresh in v0.1.

Citation identity is checked, but paper relevance, species/tissue matching and whether a statement actually supports the label still require expert review. An abstract and a URL are not ground truth. No model-generated confidence is called a probability of correctness.

Reference: [Anthropic client tool lifecycle](https://platform.claude.com/docs/en/agents-and-tools/tool-use/handle-tool-calls), [Europe PMC API](https://europepmc.org/RestfulWebService), [Ensembl lookup](https://rest.ensembl.org/documentation/info/symbol_lookup).

## Reproducibility boundaries

The manifest records input and code content hashes, runtime versions, parameters and output hashes. `uv.lock` locks the development environment. Offline notebooks rebuild numerical outputs and replay saved labels; they explicitly fail if cells or clusters change. Online reruns can retrieve different evidence and produce different model decisions. Large-scale memory/runtime behavior and cross-platform numerical equivalence are not established by the current small-data tests.

## Optional v0.2 reference methods and purification

`annotation_method` selects KNN (default), scANVI or TACCO. Optional methods consume integer raw counts over intersected genes. scANVI uses scVI followed by semi-supervised scANVI, reference/query origin batches, CPU training and saved weights. TACCO uses OT with platform iterations disabled and no recursive bisections. Full normalized assignment matrices are saved along with entropy/margin. Neither output is calibrated correctness confidence. Tissue/atlas coverage and donor-aware holdouts remain user scientific responsibilities.

`split_enabled` requires imaging single-cell input plus scANVI/TACCO. The standard SPLIT RCTD-free formula is implemented in Python; see methods.py and the packaged split-decontamination skill for exact scope. All types contribute to the denominator. Per-cell library-normalized reference counts are averaged by class. Uncovered genes retain original values. The fractional correction layer is separate from expression analysis, and raw-marker auditing remains independent of correction. Native R parity is untested. See validation-v03.md and HIERARCHY_TUTORIAL.md.

## Iterative hierarchy

See HIERARCHY_TUTORIAL.md for the complete scoring and decision specification. Hierarchical refinement follows initial annotation and precedes spatial analysis. Each round reconstructs pools from original raw counts, applies changes synchronously and clears stale subtypes. Independent atlas transfer is reference-fitted PCA/KNN in this release. Numerical stability, cycling and iteration limits are reported separately.
