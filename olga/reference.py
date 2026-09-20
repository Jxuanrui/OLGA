"""Reference-pack management: list, load and locate cell-type reference packs.

A reference pack is a directory produced by ``build_reference.py`` containing
``tau_matrix.tsv``, ``markers_strict.tsv``, ``expr_by_cluster.tsv`` and
``manifest.json``. Packs ship as derived statistics only, so raw single-cell
atlases never have to be redistributed.
"""
from __future__ import annotations

import csv
import json
import os
from pathlib import Path

_ENV = "OLGA_REFS"


def refs_root() -> Path:
    if os.environ.get(_ENV):
        return Path(os.environ[_ENV])
    bundled = Path(__file__).resolve().parent / "data" / "celltype_refs"
    if bundled.exists():
        return bundled
    # development fallback: alongside the project
    return Path(__file__).resolve().parents[1] / "data" / "reference" / "celltype_refs"


def list_packs(root: Path | None = None) -> list[dict]:
    root = root or refs_root()
    packs = []
    if not root.exists():
        return packs
    for d in sorted(root.iterdir()):
        man = d / "manifest.json"
        if man.exists():
            m = json.loads(man.read_text())
            m["_path"] = str(d)
            packs.append(m)
    return packs


def pack_names(root: Path | None = None) -> list[str]:
    return [p["name"] for p in list_packs(root)]


class Pack:
    def __init__(self, path: Path):
        self.path = Path(path)
        self.manifest = json.loads((self.path / "manifest.json").read_text())
        self.tau: dict[str, dict] = {}
        with (self.path / "tau_matrix.tsv").open() as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                self.tau[r["gene"]] = r
        self.markers: dict[str, str] = {}
        mk = self.path / "markers_strict.tsv"
        if mk.exists():
            with mk.open() as fh:
                for r in csv.DictReader(fh, delimiter="\t"):
                    self.markers[r["gene"]] = r["cell_state"]

    @property
    def name(self) -> str:
        return self.manifest["name"]

    @property
    def n_states(self) -> int:
        return int(self.manifest.get("n_cell_states", 0))

    def attribute(self, gene: str) -> dict | None:
        r = self.tau.get(gene)
        if not r:
            return None
        return {"cell_state": r["top_cell_state"],
                "tau": float(r["tau"]) if r["tau"] else None,
                "expr": float(r["top_expr"]) if r["top_expr"] else None,
                "marker_state": self.markers.get(gene)}


def load_packs(root: Path | None = None, names: list[str] | None = None) -> list[Pack]:
    root = root or refs_root()
    packs = []
    for m in list_packs(root):
        if names and m["name"] not in names:
            continue
        packs.append(Pack(Path(m["_path"])))
    return packs
