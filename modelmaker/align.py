from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable

from Bio.Align import PairwiseAligner, substitution_matrices

from .models import AlignmentResult
from .uniprot import FastaRecord, fetch_fasta_records


@dataclass(frozen=True)
class _CandidateAlignment:
    record: FastaRecord
    aligned_human: str
    aligned_zfish: str
    score: float
    identity: float
    similarity: float
    gaps: float


def _reconstruct_alignment(alignment, seq1: str, seq2: str) -> tuple[str, str]:
    a1_segments = alignment.aligned[0]
    a2_segments = alignment.aligned[1]
    out1: list[str] = []
    out2: list[str] = []
    i1 = 0
    i2 = 0
    for (s1, e1), (s2, e2) in zip(a1_segments, a2_segments):
        if i1 < s1:
            out1.append(seq1[i1:s1])
            out2.append("-" * (s1 - i1))
        if i2 < s2:
            out1.append("-" * (s2 - i2))
            out2.append(seq2[i2:s2])
        out1.append(seq1[s1:e1])
        out2.append(seq2[s2:e2])
        i1 = e1
        i2 = e2
    if i1 < len(seq1):
        out1.append(seq1[i1:])
        out2.append("-" * (len(seq1) - i1))
    if i2 < len(seq2):
        out1.append("-" * (len(seq2) - i2))
        out2.append(seq2[i2:])
    a1 = "".join(out1)
    a2 = "".join(out2)
    n = max(len(a1), len(a2))
    return a1.ljust(n, "-"), a2.ljust(n, "-")


def _score_alignment(a1: str, a2: str) -> tuple[float, float, float]:
    aligned = 0
    identical = 0
    gaps = 0
    for x, y in zip(a1, a2):
        if x == "-" or y == "-":
            gaps += 1
            continue
        aligned += 1
        if x == y:
            identical += 1
    denom = max(len(a1), 1)
    return identical / max(aligned, 1), identical / max(aligned, 1), gaps / denom


def _align_pair(seq1: str, seq2: str) -> tuple[str, str, float, float, float]:
    aligner = PairwiseAligner()
    aligner.mode = "global"
    try:
        aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
    except Exception:
        aligner.match_score = 1
        aligner.mismatch_score = -1
    aligner.open_gap_score = -10
    aligner.extend_gap_score = -0.5
    alignment = aligner.align(seq1, seq2)[0]
    a1, a2 = _reconstruct_alignment(alignment, seq1, seq2)
    identity, similarity, gaps = _score_alignment(a1, a2)
    return a1, a2, float(alignment.score), identity, similarity, gaps


def _align_candidate(h_seq: str, record: FastaRecord) -> _CandidateAlignment:
    a1, a2, score, identity, similarity, gaps = _align_pair(h_seq, record.sequence)
    return _CandidateAlignment(
        record=record,
        aligned_human=a1,
        aligned_zfish=a2,
        score=score,
        identity=identity,
        similarity=similarity,
        gaps=gaps,
    )


def _choose_candidate(h_seq: str, candidates: list[FastaRecord], preferred_gene: str = "") -> tuple[FastaRecord, str, float, float, float]:
    if not candidates:
        raise ValueError("No FASTA candidates found")
    preferred_gene = preferred_gene.strip().upper()

    if preferred_gene:
        exact = [r for r in candidates if r.gene_symbol.upper() == preferred_gene]
        if exact:
            scored = [_align_candidate(h_seq, r) for r in exact]
            best = max(scored, key=lambda x: (x.identity, x.score, -len(x.record.sequence)))
            return best.record, best.aligned_zfish, best.score, best.identity, best.gaps

    scored = [_align_candidate(h_seq, r) for r in candidates]
    best = max(scored, key=lambda x: (x.identity, x.score, -len(x.record.sequence)))
    return best.record, best.aligned_zfish, best.score, best.identity, best.gaps


def align_proteins(human_accession: str, zfish_query: str, human_gene: str = "", zfish_gene: str = "") -> AlignmentResult:
    h_records = fetch_fasta_records(human_accession)
    if not h_records:
        raise ValueError(f"No human FASTA records returned for {human_accession!r}")
    h_record = h_records[0]
    z_records = fetch_fasta_records(zfish_query)
    z_record, _, score, identity, gaps = _choose_candidate(h_record.sequence, z_records, preferred_gene=zfish_gene)

    # Need the actual alignment against the selected record.
    a1, a2, score2, identity2, similarity2, gaps2 = _align_pair(h_record.sequence, z_record.sequence)

    return AlignmentResult(
        human_gene=human_gene,
        zfish_gene=zfish_gene,
        human_accession=human_accession,
        zfish_query=zfish_query,
        human_header=h_record.header,
        zfish_header=z_record.header,
        human_sequence=h_record.sequence,
        zfish_sequence=z_record.sequence,
        aligned_human=a1,
        aligned_zfish=a2,
        score=float(score2),
        identity=identity2,
        similarity=similarity2,
        gaps=gaps2,
    )


def map_residue(aln_h: str, aln_z: str, human_position: int) -> tuple[int | None, int | None, str | None, str | None]:
    hp = 0
    zp = 0
    alignment_column = 0
    for h, z in zip(aln_h, aln_z):
        alignment_column += 1
        if h != "-":
            hp += 1
        if z != "-":
            zp += 1
        if hp == human_position and h != "-":
            return (alignment_column, zp if z != "-" else None, z if z != "-" else None, h)
    return None, None, None, None


def format_alignment(aln_h: str, aln_z: str, width: int = 60, label_h: str = "human", label_z: str = "zfish") -> str:
    lines: list[str] = []
    hp = 0
    zp = 0
    for start in range(0, len(aln_h), width):
        chunk_h = aln_h[start:start + width]
        chunk_z = aln_z[start:start + width]
        start_h = hp + 1
        start_z = zp + 1
        hp += sum(1 for c in chunk_h if c != "-")
        zp += sum(1 for c in chunk_z if c != "-")
        mid = "".join("|" if (a == b and a != "-") else " " for a, b in zip(chunk_h, chunk_z))
        lines.append(f"{label_h:<12} {start_h:>5} {chunk_h} {hp:>5}")
        lines.append(f"{'':<12} {'':>5} {mid}")
        lines.append(f"{label_z:<12} {start_z:>5} {chunk_z} {zp:>5}")
        lines.append("")
    return "\n".join(lines).rstrip()