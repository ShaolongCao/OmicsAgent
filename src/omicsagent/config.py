from pathlib import Path
from typing import Literal

from pydantic import BaseModel, ConfigDict, Field, model_validator


class RunConfig(BaseModel):
    model_config = ConfigDict(extra="forbid", allow_inf_nan=False)
    input: Path
    out: Path
    modality: Literal["scrna", "spatial"] = "scrna"
    platform: Literal["scrna", "imaging", "visium"] = "scrna"
    species: Literal["human", "mouse"] = "human"
    tissue: str = Field(min_length=1)
    counts_layer: str = "counts"  # X is only used when explicitly requested.
    seed: int = 42
    min_genes: int = Field(default=200, ge=1)
    min_cells: int = Field(default=3, ge=1)
    max_pct_mt: float = Field(default=20, ge=0, le=100)
    mt_prefix: str | None = None
    normalization: Literal["total", "area", "none"] = "total"
    area_key: str | None = None
    n_hvg: int = Field(default=2000, ge=2)
    n_pcs: int = Field(default=30, ge=2)
    n_neighbors: int = Field(default=15, ge=2)
    resolution: float = Field(default=0.5, gt=0, le=5)
    top_markers: int = Field(default=30, ge=2, le=200)
    marker_file: Path | None = None
    hierarchy: Path | None = None
    annotation_key: str | None = None
    reference: Path | None = None
    reference_label_key: str | None = None
    reference_counts_layer: str = "counts"
    min_shared_genes: int = Field(default=30, ge=3)
    annotation_method: Literal["knn", "scanvi", "tacco"] = "knn"
    scanvi_epochs: int = Field(default=50, ge=1, le=1000)
    scanvi_classifier_epochs: int = Field(default=30, ge=1, le=1000)
    scanvi_batch_size: int = Field(default=128, ge=2)
    split_enabled: bool = False
    knn_k: int = Field(default=15, ge=1)
    min_vote: float = Field(default=0.6, ge=0, le=1)
    min_marker_support: int = Field(default=1, ge=1)
    spatial_key: str = "spatial"
    sample_key: str | None = None
    coordinate_units: Literal["micron", "pixel"] | None = None
    radius: float | None = Field(default=None, gt=0)
    n_permutations: int = Field(default=99, ge=9, le=10000)
    spatial_genes: int = Field(default=10, ge=0, le=100)
    center_group_key: str | None = None
    center_groups: tuple[str, str] | None = None
    min_centers: int = Field(default=5, ge=1)
    max_edges: int = Field(default=2_000_000, ge=1)
    max_dense_elements: int = Field(default=20_000_000, ge=1000)
    online_evidence: bool = False
    model: str | None = None
    max_agent_steps: int = Field(default=20, ge=1, le=100)
    max_tool_calls: int = Field(default=50, ge=1, le=200)
    max_output_tokens: int = Field(default=2048, ge=128, le=8192)
    max_total_tokens: int = Field(default=60000, ge=1000)

    @model_validator(mode="after")
    def coherent(self) -> "RunConfig":
        if self.modality == "spatial":
            if self.platform == "scrna":
                raise ValueError("Spatial data requires platform=imaging or visium")
            if self.sample_key is None or self.coordinate_units is None or self.radius is None:
                raise ValueError("Spatial data requires sample_key, coordinate_units and radius")
        elif self.platform != "scrna":
            raise ValueError("scrna modality requires scrna platform")
        if self.normalization == "area" and self.area_key is None:
            raise ValueError("area normalization requires area_key")
        if bool(self.reference) != bool(self.reference_label_key):
            raise ValueError("reference and reference_label_key must be set together")
        if self.reference and self.annotation_key:
            raise ValueError("Choose reference transfer OR an existing annotation_key")
        if self.reference and self.annotation_method == "knn" and self.normalization != "total":
            raise ValueError("v0.1 reference transfer supports total normalization only")
        if self.annotation_method != "knn" and not self.reference:
            raise ValueError("scANVI/TACCO require a labeled raw-count reference")
        if self.split_enabled and (
            self.platform != "imaging" or self.annotation_method not in {"scanvi", "tacco"}
        ):
            raise ValueError("SPLIT core requires imaging data and scANVI/TACCO weights")
        if bool(self.center_groups) != bool(self.center_group_key):
            raise ValueError("center_group_key and center_groups must be set together")
        if self.center_groups and self.center_groups[0] == self.center_groups[1]:
            raise ValueError("Two distinct center groups are required")
        return self


def load_config(path: Path) -> RunConfig:
    """Resolve paths against the config location, not the caller's working directory."""
    config = RunConfig.model_validate_json(path.read_text())
    for name in ("input", "out", "reference", "marker_file", "hierarchy"):
        value = getattr(config, name)
        if value is not None:
            setattr(config, name, (path.resolve().parent / value).resolve())
    return config
