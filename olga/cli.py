"""Command-line interface for OLGA."""
from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

from . import __version__
from .attribution import attribute_gene
from .buildref import build_reference
from .outputs import plot_evidence_tiers, plot_gene_celltype, write_chains, write_manifest
from .reference import list_packs, load_packs
from .verify import verify as run_verify


def cmd_list_reference(args) -> int:
    packs = list_packs(Path(args.refs) if args.refs else None)
    if not packs:
        print("no reference packs found; run `olga init` or set OLGA_REFS")
        return 1
    print(f"{'name':28s} {'states':>7s} {'genes':>7s}")
    for m in packs:
        print(f"{m['name']:28s} {m.get('n_cell_states', '?'):>7} {m.get('n_genes', '?'):>7}")
    return 0


def _run_attribution(records, refs, names, out) -> int:
    """records: list of (gene, meta) — meta carries optional trait/locus context."""
    packs = load_packs(Path(refs) if refs else None, names)
    if not packs:
        print("error: no reference packs available", file=sys.stderr)
        return 2
    out = Path(out)
    results = []
    for gene, meta in records:
        r = attribute_gene(gene, packs)
        r.update(meta)
        results.append(r)
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
                    "n_genes": len(records), "evidence_tiers": tiers,
                    "outputs": ["chains.tsv", "celltypes/*.celltype.{png,pdf}",
                                "evidence_tiers.png"]},
                   out / "run_manifest.json")
    print(f"attributed {len(records)} genes with {len(packs)} reference packs -> {out}")
    for r in results:
        ctx = f" [{r['trait_id']} {r['locus_id']}]" if r.get("trait_id") else ""
        print(f"  {r['gene']:16s}{ctx} {str(r['consensus_state']):32s} {r['evidence_tier']}")
    return 0


def cmd_attribute(args) -> int:
    genes = [g.strip() for g in args.genes.split(",") if g.strip()]
    return _run_attribution([(g, {}) for g in genes], args.refs, args.reference, args.out)


def read_effector_records(path: str) -> tuple[list, int]:
    """Parse an effector-gene TSV. Whole-chain mode when the header carries
    trait_id/locus_id columns; returns (records, n_rows_skipped_no_gene)."""
    header = open(path).readline().rstrip("\n").split("\t")
    records, skipped = [], 0
    if "effector_gene" in header:
        with open(path) as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                gene = (r.get("effector_gene") or "").strip()
                if not gene or gene.lower() == "na":
                    skipped += 1
                    continue
                meta = {k: r[k] for k in ("trait_id", "locus_id", "lead_snp")
                        if k in header}
                records.append((gene, meta))
    elif "trait_id" in header:
        raise SystemExit("genes-file looks like a trait table but has no "
                         "effector_gene column")
    else:
        with open(path) as fh:
            for line in fh:
                line = line.strip()
                if not line or line.startswith("#"):
                    continue
                records.append((line.split("\t")[0].strip(), {}))
    return records, skipped


def cmd_run(args) -> int:
    records, skipped = ([], 0)
    if args.genes_file:
        records, skipped = read_effector_records(args.genes_file)
    if args.genes:
        records += [(g.strip(), {}) for g in args.genes.split(",") if g.strip()]
    if not records:
        print("error: provide effector genes via --genes or --genes-file", file=sys.stderr)
        return 2
    if skipped:
        print(f"note: {skipped} table rows without an effector gene skipped",
              file=sys.stderr)
    return _run_attribution(records, args.refs, args.reference, args.out)


def cmd_verify(args) -> int:
    return run_verify()


def cmd_build_reference(args) -> int:
    labels = args.label or []
    atlases = list(zip(args.atlas, labels + [None] * (len(args.atlas) - len(labels))))
    m = build_reference(atlases, args.name, args.out, args.cluster_col)
    print(f"{m['name']}: cells={m['n_cells']} genes={m['n_genes']} "
          f"states={m['n_cell_states']} cluster_col={m['cluster_column']} -> {args.out}")
    return 0


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
                    help="TSV of effector genes (an OLGA effector table with "
                         "trait_id/locus_id/effector_gene columns carries the "
                         "whole chain through to chains.tsv; otherwise the "
                         "first column is gene symbols)")
    rn.add_argument("--reference", nargs="*", default=None)
    rn.add_argument("--refs", default=None)
    rn.add_argument("--out", required=True)
    rn.set_defaults(func=cmd_run)

    vf = sub.add_parser("verify", help="self-check packs, lineage rules and golden attributions")
    vf.set_defaults(func=cmd_verify)

    br = sub.add_parser("build-reference", help="build a cell-type reference pack from annotated h5ad atlas(es)")
    br.add_argument("--atlas", action="append", required=True,
                    help="input .h5ad (repeat to merge atlases on the gene union)")
    br.add_argument("--label", action="append", default=None,
                    help="state-name prefix per atlas, aligned order (e.g. --label Epi --label LP)")
    br.add_argument("--name", required=True, help="pack name")
    br.add_argument("--out", required=True, help="output pack directory")
    br.add_argument("--cluster-col", default=None,
                    help="obs column with cell states (auto-detected if omitted)")
    br.set_defaults(func=cmd_build_reference)
    return p


def main(argv=None) -> int:
    args = build_parser().parse_args(argv)
    return args.func(args)


if __name__ == "__main__":
    sys.exit(main())
