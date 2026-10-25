import argparse
import json
import os
from pathlib import Path
from typing import Any


def main() -> None:
    parser = argparse.ArgumentParser(
        description="OmicsAgent v0.3 — scRNA-seq and spatial workflows"
    )
    sub = parser.add_subparsers(dest="command", required=True)
    for name in ("run", "agent"):
        command = sub.add_parser(name)
        command.add_argument("--config", type=Path, required=True)
    demo = sub.add_parser("demo", help="Generate synthetic integration data and example configs")
    demo.add_argument("--out", type=Path, required=True)
    hdemo = sub.add_parser("hierarchy-demo", help="Generate independent synthetic atlas and iterative annotation examples")
    hdemo.add_argument("--out", type=Path, required=True)
    inspect = sub.add_parser("inspect", help="Read-only h5ad schema inspection")
    inspect.add_argument("path", type=Path)
    sub.add_parser(
        "doctor", help="Check runtime and presence of API configuration; never prints keys"
    )
    fetch = sub.add_parser("fetch-pbmc3k", help="Download the public Scanpy raw PBMC3k dataset")
    fetch.add_argument("--out", type=Path, required=True)
    ingest = sub.add_parser("import-10x", help="Convert 10x H5 plus barcode-indexed CSV to h5ad")
    ingest.add_argument("--matrix", type=Path, required=True)
    ingest.add_argument("--metadata", type=Path)
    ingest.add_argument("--out", type=Path, required=True)
    ingest.add_argument("--id-column", default="cell_id")
    ingest.add_argument("--x-column")
    ingest.add_argument("--y-column")
    args = parser.parse_args()
    result: dict[str, Any]
    try:
        if args.command == "hierarchy-demo":
            from .demo import make_hierarchy_demo

            make_hierarchy_demo(args.out)
            result = {"status": "created", "out": str(args.out)}
        elif args.command == "demo":
            from .demo import make_demo

            make_demo(args.out)
            result = {"status": "created", "out": str(args.out)}
        elif args.command == "import-10x":
            from .ingest import import_10x

            result = import_10x(
                args.matrix, args.metadata, args.out, args.id_column, args.x_column, args.y_column
            )
        elif args.command == "inspect":
            import anndata as ad

            obj = ad.read_h5ad(args.path, backed="r")
            try:
                result = {
                    "shape": list(obj.shape),
                    "obs_columns": list(obj.obs.columns),
                    "layers": list(obj.layers),
                    "obsm": list(obj.obsm),
                    "genes_example": list(obj.var_names[:10]),
                }
            finally:
                obj.file.close()
        elif args.command == "doctor":
            import importlib.metadata
            import platform

            optional = {}
            for package in ("scvi-tools", "tacco", "torch", "harmonypy"):
                try:
                    optional[package] = importlib.metadata.version(package)
                except importlib.metadata.PackageNotFoundError:
                    optional[package] = "not installed"
            result = {
                "optional_dependencies": optional,
                "split_backend": "Python standard core; not native R SPLIT",
                "python": platform.python_version(),
                "dependencies": {
                    p: importlib.metadata.version(p) for p in ["scanpy", "anndata", "anthropic"]
                },
                "anthropic_key_present": bool(os.environ.get("ANTHROPIC_API_KEY")),
            }
        elif args.command == "fetch-pbmc3k":
            import scanpy as sc

            if args.out.exists():
                raise ValueError("Refusing to overwrite existing input")
            args.out.parent.mkdir(parents=True, exist_ok=True)
            sc.settings.datasetdir = args.out.parent / "scanpy-cache"
            obj = sc.datasets.pbmc3k()
            obj.layers["counts"] = obj.X.copy()
            obj.write_h5ad(args.out)
            result = {"status": "downloaded", "shape": list(obj.shape), "out": str(args.out)}
        else:
            from .config import load_config

            config = load_config(args.config)
            if args.command == "agent":
                from .agent import run_agent

                result = run_agent(config)
            else:
                from .workflow import run_pipeline

                result = run_pipeline(config)
        print(json.dumps(result, indent=2))
        if result.get("status") == "incomplete":
            raise SystemExit(2)
    except (ValueError, RuntimeError, OSError) as exc:
        parser.exit(1, f"OmicsAgent: {exc}\n")
