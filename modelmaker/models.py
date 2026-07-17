from __future__ import annotations

from dataclasses import asdict, dataclass, field
from typing import Any, Dict, List, Optional


@dataclass(frozen=True)
class MutationSpec:
    gene: str
    mutation: str
    position: int
    numbering_offset: int = 0

    @property
    def canonical_position(self) -> int:
        return self.position + self.numbering_offset

    @property
    def ref_aa(self) -> str:
        return self.mutation[0]

    @property
    def alt_aa(self) -> str:
        return self.mutation[-1]


@dataclass(frozen=True)
class TargetSpec:
    gene: str
    human: str
    zfish: str


@dataclass(frozen=True)
class OrthologHit:
    gene: str
    source: str
    header: str
    sequence: str
    url: str


@dataclass(frozen=True)
class AlignmentResult:
    human_gene: str
    zfish_gene: str
    human_accession: str
    zfish_query: str
    human_header: str
    zfish_header: str
    human_sequence: str
    zfish_sequence: str
    aligned_human: str
    aligned_zfish: str
    score: float
    identity: float
    similarity: float
    gaps: float


@dataclass(frozen=True)
class ProteinInfo:
    """Core protein-level object for ModelMaker.

    ProteinInfo is intended to represent information intrinsic to a protein
    sequence, while remaining independent from genome-editing or variant design
    logic.
    """

    protein_id: str
    gene_symbol: Optional[str] = None
    species: Optional[str] = None
    sequence: Optional[str] = None

    # Alignment / similarity metadata
    best_match_id: Optional[str] = None
    percent_identity: Optional[float] = None
    alignment_coverage: Optional[float] = None
    evalue: Optional[float] = None
    similarity: Optional[float] = None
    gaps: Optional[float] = None
    alignment_score: Optional[float] = None
    aligned_sequence: Optional[str] = None

    # Basic biology annotations
    domains: List[str] = field(default_factory=list)
    motifs: List[str] = field(default_factory=list)
    biological_role: Optional[str] = None

    # Provenance / bookkeeping
    source: Optional[str] = None
    notes: Optional[str] = None

    def to_dict(self) -> Dict[str, Any]:
        return asdict(self)


@dataclass(frozen=True)
class MappingResult:
    mutation: MutationSpec
    human_canonical_position: int
    human_residue: str
    zfish_position: Optional[int]
    zfish_residue: Optional[str]
    conserved: bool
    human_accession: str
    zfish_query: str


@dataclass(frozen=True)
class NearbyFunctionalSite:
    found: bool
    within_window: bool
    human_position: Optional[int]
    zfish_position: Optional[int]
    alignment_column: Optional[int]
    human_residue: Optional[str]
    zfish_residue: Optional[str]
    distance: Optional[int]
    side: str = ""
    note: str = ""


@dataclass(frozen=True)
class EditCandidate:
    editor: str
    spacer: str
    pam: str
    strand: str
    target_index: int
    window_start: int
    window_end: int
    bystanders: List[Dict[str, Any]]
    notes: str = ""


@dataclass(frozen=True)
class EditabilityResult:
    editable: bool
    recommended_strategy: str
    confidence: str
    reason_codes: List[str]
    pam_present: bool
    pam_sequence: str
    editing_window_covers_target: bool
    bystander_risk: str
    codon_reachability: Dict[str, bool]
    sgRNA_candidate: Optional[EditCandidate] = None
    ssODN_candidate: Optional[Dict[str, Any]] = None
    pegRNA_candidate: Optional[Dict[str, Any]] = None


@dataclass(frozen=True)
class ResultBundle:
    alignment: Optional[AlignmentResult]
    mapping: Optional[MappingResult]
    editability: Optional[EditabilityResult]

    def to_dict(self) -> Dict[str, Any]:
        return asdict(self)
