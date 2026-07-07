from __future__ import annotations

import argparse
import re
from pathlib import Path

from .align import align_proteins, format_alignment, map_residue
from .editability import assess_editability
from .functional import find_nearby_functional_site
from .io import read_mutations_csv, read_targets_csv, write_csv, write_json


def _header_gene(header: str) -> str:
    m = re.search(r"\bGN=([^\s]+)", header)
    return m.group(1) if m else ""


def _write_alignment_outputs(outdir: Path, gene: str, aln, width: int) -> None:
    outdir.mkdir(parents=True, exist_ok=True)
    (outdir / f"{gene}_alignment.txt").write_text(
        format_alignment(aln.aligned_human, aln.aligned_zfish, width=width, label_h="human", label_z="zfish") + "\n",
        encoding="utf-8",
    )
    write_json(outdir / f"{gene}_alignment.json", aln)


def _map_one_target(target_gene: str, human: str, zfish: str, variants, outdir: Path, width: int, nearby_window: int):
    aln = align_proteins(human, zfish, human_gene=target_gene, zfish_gene=target_gene)
    target_dir = outdir / target_gene.lower()
    _write_alignment_outputs(target_dir, target_gene, aln, width)

    rows: list[dict[str, object]] = []
    summary_lines: list[str] = []
    human_header_gene = _header_gene(aln.human_header)
    zfish_header_gene = _header_gene(aln.zfish_header)

    for mut in variants:
        if mut.gene.upper() != target_gene.upper():
            continue
        mature_pos = mut.position
        canonical_pos = mut.canonical_position
        aln_col, zpos, zaa, haa = map_residue(aln.aligned_human, aln.aligned_zfish, canonical_pos)
        conserved = bool(zaa and zaa == mut.ref_aa)
        nearby = None
        if not conserved:
            nearby = find_nearby_functional_site(aln.aligned_human, aln.aligned_zfish, canonical_pos, window=nearby_window)

        row = {
            "gene": mut.gene,
            "mutation": mut.mutation,
            "reported_position": mature_pos,
            "numbering_offset": mut.numbering_offset,
            "human_mature_position": mature_pos,
            "human_canonical_position": canonical_pos,
            "alignment_column": aln_col if aln_col is not None else "",
            "human_residue": haa or mut.ref_aa,
            "zfish_position": zpos if zpos is not None else "",
            "zfish_residue": zaa if zaa is not None else "",
            "human_header_gene": human_header_gene,
            "zfish_header_gene": zfish_header_gene,
            "conserved": "yes" if conserved else "no",
            "nearby_functional_site_found": "yes" if nearby and nearby.found else "no",
            "nearby_functional_site_within_window": "yes" if nearby and nearby.within_window else "no",
            "nearby_functional_site_human_position": nearby.human_position if nearby and nearby.found else "",
            "nearby_functional_site_zfish_position": nearby.zfish_position if nearby and nearby.found else "",
            "nearby_functional_site_alignment_column": nearby.alignment_column if nearby and nearby.found else "",
            "nearby_functional_site_human_residue": nearby.human_residue if nearby and nearby.found else "",
            "nearby_functional_site_zfish_residue": nearby.zfish_residue if nearby and nearby.found else "",
            "nearby_functional_site_distance": nearby.distance if nearby and nearby.found else "",
            "nearby_functional_site_side": nearby.side if nearby and nearby.found else "",
            "nearby_functional_site_note": nearby.note if nearby and nearby.found else ("Exact residue conserved; no nearby-site search needed." if conserved else ""),
        }
        rows.append(row)
        summary = (
            f"{mut.gene:<7} {mut.mutation:<8} mature {mature_pos} canonical {canonical_pos} "
            f"align_col {aln_col} human {haa} -> zfish {zpos} {zaa} conserved={'yes' if conserved else 'no'}"
        )
        if nearby and nearby.found:
            summary += f" | nearby {nearby.human_position} {nearby.human_residue} -> {nearby.zfish_position} {nearby.zfish_residue}"
        summary_lines.append(summary)
    return aln, rows, summary_lines


def cmd_map(args: argparse.Namespace) -> int:
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    if args.targets and args.variants:
        targets = read_targets_csv(Path(args.targets))
        variants = read_mutations_csv(Path(args.variants))
        all_rows: list[dict[str, object]] = []
        all_summary: list[str] = []
        for target in targets:
            _, rows, summary_lines = _map_one_target(
                target.gene,
                target.human,
                target.zfish,
                variants,
                outdir,
                args.width,
                args.nearby_window,
            )
            all_rows.extend(rows)
            all_summary.extend([f"[{target.gene}] {line}" for line in summary_lines])
        write_csv(outdir / "variant_map.csv", all_rows)
        (outdir / "summary.txt").write_text("\n".join(all_summary) + "\n", encoding="utf-8")
        return 0

    if not args.human or not args.zfish or not args.mutations:
        raise SystemExit("Provide either --targets/--variants or --human/--zfish/--mutations")

    aln, rows, summary_lines = _map_one_target(
        args.human_gene or args.zfish_gene or _header_gene(args.human) or "GENE",
        args.human,
        args.zfish,
        read_mutations_csv(Path(args.mutations)),
        outdir,
        args.width,
        args.nearby_window,
    )
    write_csv(outdir / "variant_map.csv", rows)
    (outdir / "summary.txt").write_text("\n".join(summary_lines) + "\n", encoding="utf-8")
    write_json(outdir / "alignment.json", aln)
    (outdir / "alignment.txt").write_text(
        format_alignment(aln.aligned_human, aln.aligned_zfish, width=args.width, label_h="human", label_z="zfish") + "\n",
        encoding="utf-8",
    )
    return 0


def cmd_design_base(args: argparse.Namespace) -> int:
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    cds = Path(args.cds).read_text(encoding="utf-8")
    cds = "".join([c for c in cds.splitlines() if not c.startswith(">")]).replace(" ", "").replace("\t", "").upper()
    result = assess_editability(cds, args.target_index, args.ref_base, args.alt_base)
    write_json(outdir / "editability.json", result)
    if result.sgRNA_candidate:
        write_json(outdir / "sgRNA_candidate.json", result.sgRNA_candidate)
    return 0


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(prog="varianttool")
    sub = p.add_subparsers(dest="cmd", required=True)

    m = sub.add_parser("map", help="Download proteins, align orthologs, and map residues")
    m.add_argument("--human", default="", help="Human UniProt accession")
    m.add_argument("--zfish", default="", help="Zebrafish UniProt accession or query")
    m.add_argument("--mutations", default="", help="CSV with gene,mutation,position,numbering_offset")
    m.add_argument("--targets", default="", help="CSV with gene,human,zfish for multi-target runs")
    m.add_argument("--variants", default="", help="CSV with gene,mutation,position,numbering_offset for multi-target runs")
    m.add_argument("--outdir", default="out", help="Output folder")
    m.add_argument("--width", type=int, default=60)
    m.add_argument("--nearby-window", type=int, default=12, help="Search radius in amino acids for nearby conserved sites")
    m.add_argument("--human-gene", default="")
    m.add_argument("--zfish-gene", default="")
    m.set_defaults(func=cmd_map)

    d = sub.add_parser("design-base", help="Assess base-editability and emit a candidate sgRNA if possible")
    d.add_argument("--cds", required=True, help="Coding DNA sequence in FASTA/plaintext form")
    d.add_argument("--target-index", type=int, required=True, help="0-based nucleotide index in CDS")
    d.add_argument("--ref-base", required=True)
    d.add_argument("--alt-base", required=True)
    d.add_argument("--outdir", default="out")
    d.set_defaults(func=cmd_design_base)

    return p


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    return args.func(args)


if __name__ == "__main__":
    raise SystemExit(main())
