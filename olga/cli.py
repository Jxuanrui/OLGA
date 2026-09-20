"""Command-line interface for OLGA."""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

from . import __version__
from .attribution import attribute_gene
from .outputs import plot_evidence_tiers, plot_gene_celltype, write_chains, write_manifest
from .reference import list_packs, load_packs


def cmd_list_reference(args) -> int:
    packs = list_packs(Path(args.refs) if args.refs else None)
    if not packs:
        print("no reference packs found; run `olga init` or set OLGA_REFS")
        return 1
    print(f"{'name':28s} {'states':>7s} {'genes':>7s}")
    for m in packs:
        print(f"{m['name']:28s} {m.get('n_cell_states', '?'):>7} {m.get('n_genes', '?'):>7}")
    return 0


def _run_attribution(genes, refs, names, out) -> int:
    packs = load_packs(Path(refs) if refs else None, names)
    if not packs:
        print("error: no reference packs available", file=sys.stderr)
        return 2
    out = Path(out)
    results = [attribute_gene(g, packs) for g in genes]
    write_chains(results, out / "chains.tsv")
    (out / "celltypes").mkdir(parents=True, exist_ok=True)
    for r in results:
        if r["n_refs_covering"]:
            plot_gene_celltype(r, out / "celltypes" / f"{r['gene']}.celltype.png",
                               out / "celltypes" / f"{r['gene']}.celltype.pdf")
    plot_evidence_tiers(results, out / "evidence_tiers.png")
    tiers = {}
    for r in results:
        tiers[r["evidence_tier"]] = tiers.get(r["evidence_tier"], 0) + 1
    write_manifest({"tool": "olga", "version": __version__,
                    "reference_packs": [p.name for p in packs],
                    "n_genes": len(genes), "evidence_tiers": tiers,
                    "outputs": ["chains.tsv", "celltypes/*.celltype.{png,pdf}",
                                "evidence_tiers.png"]},
                   out / "run_manifest.json")
    print(f"attributed {len(genes)} genes with {len(packs)} reference packs -> {out}")
    for r in results:
        print(f"  {r['gene']:12s} {str(r['consensus_state']):32s} {r['evidence_tier']}")
    return 0


def cmd_attribute(args) -> int:
    genes = [g.strip() for g in args.genes.split(",") if g.strip()]
    return _run_attribution(genes, args.refs, args.reference, args.out)


def cmd_run(args) -> int:
    genes = []
    if args.genes_file:
        with open(args.genes_file) as fh:
            for line in fh:
                line = line.strip()
                if not line or line.startswith("#"):
                    continue
                genes.append(line.split("\t")[0].strip())
    if args.genes:
        genes += [g.strip() for g in args.genes.split(",") if g.strip()]
    if not genes:
        print("error: provide effector genes via --genes or --genes-file", file=sys.stderr)
        return 2
    return _run_attribution(genes, args.refs, args.reference, args.out)


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(prog="olga", description=__doc__)
    p.add_argument("--version", action="version", version=f"olga {__version__}")
    sub = p.add_subparsers(dest="cmd", required=True)

    lp = sub.add_parser("list-reference", help="list bundled cell-type reference packs")
    lp.add_argument("--refs", default=None)
    lp.set_defaults(func=cmd_list_reference)

    at = sub.add_parser("attribute", help="attribute effector genes to cell types")
    at.add_argument("--genes", required=True, help="comma-separated gene symbols")
    at.add_argument("--reference", nargs="*", default=None, help="subset of pack names")
    at.add_argument("--refs", default=None, help="reference root directory")
    at.add_argument("--out", required=True)
    at.set_defaults(func=cmd_attribute)

    rn = sub.add_parser("run", help="run the attribution pipeline on effector genes")
    rn.add_argument("--genes", default=None, help="comma-separated gene symbols")
    rn.add_argument("--genes-file", default=None,
                    help="TSV whose first column is gene symbols")
    rn.add_argument("--reference", nargs="*", default=None)
    rn.add_argument("--refs", default=None)
    rn.add_argument("--out", required=True)
    rn.set_defaults(func=cmd_run)
    return p


def main(argv=None) -> int:
    args = build_parser().parse_args(argv)
    return args.func(args)


if __name__ == "__main__":
    sys.exit(main())
