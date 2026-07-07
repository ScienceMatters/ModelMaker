from __future__ import annotations

from dataclasses import asdict
from typing import Dict, Iterable, List, Optional

from .models import NearbyFunctionalSite


def iter_aligned_positions(aln_h: str, aln_z: str) -> List[Dict[str, object]]:
    """Return a list of aligned positions with residue numbering in both species."""
    hp = 0
    zp = 0
    out: List[Dict[str, object]] = []
    for col, (h, z) in enumerate(zip(aln_h, aln_z), start=1):
        if h != '-':
            hp += 1
        if z != '-':
            zp += 1
        if h == '-' or z == '-':
            conserved = False
        else:
            conserved = h == z
        out.append(
            {
                'alignment_column': col,
                'human_position': hp if h != '-' else None,
                'zfish_position': zp if z != '-' else None,
                'human_residue': None if h == '-' else h,
                'zfish_residue': None if z == '-' else z,
                'conserved': conserved,
            }
        )
    return out


def map_alignment_position(aln_h: str, aln_z: str, human_position: int) -> Dict[str, object]:
    hp = 0
    zp = 0
    for col, (h, z) in enumerate(zip(aln_h, aln_z), start=1):
        if h != '-':
            hp += 1
        if z != '-':
            zp += 1
        if hp == human_position and h != '-':
            return {
                'alignment_column': col,
                'zfish_position': zp if z != '-' else None,
                'zfish_residue': None if z == '-' else z,
                'human_residue': h,
            }
    return {
        'alignment_column': None,
        'zfish_position': None,
        'zfish_residue': None,
        'human_residue': None,
    }


def find_nearby_functional_site(
    aln_h: str,
    aln_z: str,
    human_position: int,
    window: int = 12,
) -> NearbyFunctionalSite:
    """Find the nearest conserved aligned residue around a target human position.

    This is a practical proxy for a nearby functional site when the exact residue
    is not conserved. It searches within +/- window human residues first, then
    falls back to the nearest conserved aligned residue anywhere in the protein.
    """
    aligned = iter_aligned_positions(aln_h, aln_z)
    candidates: List[Dict[str, object]] = [
        p for p in aligned
        if p['conserved'] and p['human_position'] is not None and p['human_position'] != human_position
    ]
    if not candidates:
        return NearbyFunctionalSite(
            found=False,
            within_window=False,
            human_position=None,
            zfish_position=None,
            alignment_column=None,
            human_residue=None,
            zfish_residue=None,
            distance=None,
            note='No conserved aligned residues found in the protein.',
        )

    within = [p for p in candidates if abs(int(p['human_position']) - human_position) <= window]
    pool = within if within else candidates
    best = min(
        pool,
        key=lambda p: (
            abs(int(p['human_position']) - human_position),
            0 if int(p['human_position']) < human_position else 1,
            int(p['alignment_column']),
        ),
    )

    distance = abs(int(best['human_position']) - human_position)
    side = 'upstream' if int(best['human_position']) < human_position else 'downstream'
    note = (
        f"Nearest conserved residue {side} of target at human {best['human_position']} / "
        f"zebrafish {best['zfish_position']} (alignment column {best['alignment_column']})."
    )
    if not within:
        note += f" Outside the +/-{window} aa search window; reported as the best conserved fallback."

    return NearbyFunctionalSite(
        found=True,
        within_window=bool(within),
        human_position=int(best['human_position']),
        zfish_position=int(best['zfish_position']) if best['zfish_position'] is not None else None,
        alignment_column=int(best['alignment_column']),
        human_residue=str(best['human_residue']) if best['human_residue'] is not None else None,
        zfish_residue=str(best['zfish_residue']) if best['zfish_residue'] is not None else None,
        distance=distance,
        side=side,
        note=note,
    )


def nearby_site_to_dict(site: NearbyFunctionalSite) -> Dict[str, object]:
    return asdict(site)
