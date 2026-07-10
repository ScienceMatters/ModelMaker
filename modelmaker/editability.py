from __future__ import annotations

from dataclasses import asdict
from typing import Dict, List, Optional, Tuple

from Bio.Seq import Seq

from .models import EditCandidate, EditabilityResult

# Standard genetic code.
CODON_TABLE = {
    'TTT':'F','TTC':'F','TTA':'L','TTG':'L',
    'CTT':'L','CTC':'L','CTA':'L','CTG':'L',
    'ATT':'I','ATC':'I','ATA':'I','ATG':'M',
    'GTT':'V','GTC':'V','GTA':'V','GTG':'V',
    'TCT':'S','TCC':'S','TCA':'S','TCG':'S',
    'CCT':'P','CCC':'P','CCA':'P','CCG':'P',
    'ACT':'T','ACC':'T','ACA':'T','ACG':'T',
    'GCT':'A','GCC':'A','GCA':'A','GCG':'A',
    'TAT':'Y','TAC':'Y','TAA':'*','TAG':'*',
    'CAT':'H','CAC':'H','CAA':'Q','CAG':'Q',
    'AAT':'N','AAC':'N','AAA':'K','AAG':'K',
    'GAT':'D','GAC':'D','GAA':'E','GAG':'E',
    'TGT':'C','TGC':'C','TGA':'*','TGG':'W',
    'CGT':'R','CGC':'R','CGA':'R','CGG':'R',
    'AGT':'S','AGC':'S','AGA':'R','AGG':'R',
    'GGT':'G','GGC':'G','GGA':'G','GGG':'G',
}
AA_TO_CODONS: Dict[str, List[str]] = {}
for codon, aa in CODON_TABLE.items():
    AA_TO_CODONS.setdefault(aa, []).append(codon)

TRANSITION = {('A', 'G'), ('G', 'A'), ('C', 'T'), ('T', 'C')}
PURINES = {'A', 'G'}
PYRIMIDINES = {'C', 'T'}


def reverse_complement(seq: str) -> str:
    return str(Seq(seq).reverse_complement())


def _possible_single_nt_changes(ref_codon: str, alt_aa: str) -> List[tuple[int, str, str, str]]:
    hits = []
    ref_codon = ref_codon.upper()
    for alt_codon in AA_TO_CODONS.get(alt_aa, []):
        diffs = [(i, ref_codon[i], alt_codon[i]) for i in range(3) if ref_codon[i] != alt_codon[i]]
        if len(diffs) == 1:
            i, r, a = diffs[0]
            hits.append((i, r, a, alt_codon))
    return hits


def codon_editability(ref_codon: str, alt_aa: str) -> Dict[str, bool]:
    ref_codon = ref_codon.upper()
    abe = False
    cbe = False
    prime = False
    for _, r, a, _ in _possible_single_nt_changes(ref_codon, alt_aa):
        if (r, a) in {('A', 'G'), ('T', 'C')}:
            abe = True
        if (r, a) in {('C', 'T'), ('G', 'A')}:
            cbe = True
        prime = True
    return {"abe_possible": abe, "cbe_possible": cbe, "prime_possible": prime}


def infer_bystanders(window_seq: str, target_offset: int, editor: str) -> List[Dict[str, object]]:
    window_seq = window_seq.upper()
    out = []
    if editor.startswith("ABE"):
        editable = {'A', 'T'}
    elif editor.startswith("CBE"):
        editable = {'C', 'G'}
    else:
        editable = set()
    for i, base in enumerate(window_seq):
        if i == target_offset:
            continue
        if base in editable:
            out.append({"offset": i, "base": base, "risk": "possible bystander"})
    return out


def find_candidate_guides(
    cds: str,
    target_index: int,
    ref_base: str,
    alt_base: str,
    editor: str = "ABE",
    pam: str = "NGG",
    window: tuple[int, int] = (4, 8),
) -> List[EditCandidate]:
    """Find candidate spacers for base editing around a target base.

    cds is expected to be 5'->3' coding sequence.
    target_index is 0-based in cds coordinates.
    """
    cds = cds.upper().replace("\n", "").replace(" ", "")
    ref_base = ref_base.upper()
    alt_base = alt_base.upper()
    candidates: List[EditCandidate] = []
    protospacer_len = 20

    if target_index < 0 or target_index >= len(cds):
        return candidates

    # plus strand search: PAM downstream of spacer
    for spacer_start in range(max(0, target_index - window[1]), min(target_index - window[0] + 1, len(cds) - protospacer_len + 1)):
        spacer = cds[spacer_start:spacer_start + protospacer_len]
        pam_start = spacer_start + protospacer_len
        pam_seq = cds[pam_start:pam_start + 3]
        if len(pam_seq) < 3:
            continue
        # very simple NGG or NNN checking
        if pam == "NGG" and not (pam_seq[1:] == "GG"):
            continue
        if target_index < spacer_start or target_index >= spacer_start + protospacer_len:
            continue
        offset = target_index - spacer_start
        if not (window[0] <= offset <= window[1]):
            continue
        window_seq = spacer[max(0, offset - 4): min(len(spacer), offset + 5)]
        bystanders = infer_bystanders(window_seq, min(4, offset), editor)
        candidates.append(EditCandidate(
            editor=editor,
            spacer=spacer,
            pam=pam_seq,
            strand="+",
            target_index=target_index,
            window_start=spacer_start + window[0],
            window_end=spacer_start + window[1],
            bystanders=bystanders,
            notes="plus-strand candidate",
        ))

    # minus strand search: PAM upstream on genomic plus-strand coordinates
    rc = reverse_complement(cds)
    rc_target = len(cds) - 1 - target_index
    for spacer_start in range(max(0, rc_target - window[1]), min(rc_target - window[0] + 1, len(rc) - protospacer_len + 1)):
        spacer = rc[spacer_start:spacer_start + protospacer_len]
        pam_start = spacer_start + protospacer_len
        pam_seq = rc[pam_start:pam_start + 3]
        if len(pam_seq) < 3:
            continue
        if pam == "NGG" and not (pam_seq[1:] == "GG"):
            continue
        if rc_target < spacer_start or rc_target >= spacer_start + protospacer_len:
            continue
        offset = rc_target - spacer_start
        if not (window[0] <= offset <= window[1]):
            continue
        window_seq = spacer[max(0, offset - 4): min(len(spacer), offset + 5)]
        bystanders = infer_bystanders(window_seq, min(4, offset), editor)
        candidates.append(EditCandidate(
            editor=editor,
            spacer=spacer,
            pam=pam_seq,
            strand="-",
            target_index=target_index,
            window_start=spacer_start + window[0],
            window_end=spacer_start + window[1],
            bystanders=bystanders,
            notes="minus-strand candidate",
        ))

    return candidates


def assess_editability(
    cds: str,
    target_index: int,
    ref_base: str,
    alt_base: str,
    codon_start: Optional[int] = None,
) -> EditabilityResult:
    cds = cds.upper().replace("\n", "").replace(" ", "")
    ref_base = ref_base.upper()
    alt_base = alt_base.upper()

    if codon_start is None:
        codon_start = (target_index // 3) * 3
    ref_codon = cds[codon_start: codon_start + 3]
    reach = codon_editability(ref_codon, alt_base)

    # Basic editor classification.
    if (ref_base, alt_base) in {('A', 'G'), ('T', 'C')}:
        strategy = "ABE"
    elif (ref_base, alt_base) in {('C', 'T'), ('G', 'A')}:
        strategy = "CBE"
    else:
        strategy = "PrimeEditing" if reach["prime_possible"] else "None"

    candidates = find_candidate_guides(cds, target_index, ref_base, alt_base, editor=strategy if strategy in {"ABE", "CBE"} else "ABE")
    sg = candidates[0] if candidates else None

    pam_present = bool(candidates)
    editable = strategy in {"ABE", "CBE", "PrimeEditing"}
    confidence = "high" if sg else ("medium" if strategy in {"ABE", "CBE"} else "low")
    bystander_risk = "high" if sg and len(sg.bystanders) >= 2 else ("medium" if sg and len(sg.bystanders) == 1 else "low")
    reason_codes = []
    if not pam_present:
        reason_codes.append("NO_PAM_IN_EDITING_WINDOW")
    if not editable:
        reason_codes.append("NO_SINGLE_BASE_PATH")

    ssodn = None
    if strategy == "None":
        ssodn = {
            "enabled": False,
            "reason": "No base-edit or prime-edit path inferred from this simple codon check.",
        }

    return EditabilityResult(
        editable=editable,
        recommended_strategy=strategy,
        confidence=confidence,
        reason_codes=reason_codes,
        pam_present=pam_present,
        pam_sequence=(sg.pam if sg else ""),
        editing_window_covers_target=bool(sg),
        bystander_risk=bystander_risk,
        codon_reachability=reach,
        sgRNA_candidate=sg,
        ssODN_candidate=ssodn,
        pegRNA_candidate=None,
    )
