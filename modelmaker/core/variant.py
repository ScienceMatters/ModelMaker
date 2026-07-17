"""
Core biological data model for ModelMaker.

The Variant is the central object that flows through the entire
ModelMaker pipeline. Modules enrich the Variant with biological
knowledge but never perform calculations inside this file.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from datetime import datetime, UTC
from uuid import UUID, uuid4


# ---------------------------------------------------------------------
# Base classes
# ---------------------------------------------------------------------
@dataclass
class EvidenceContainer:
    """
    Base class for biological knowledge.
    """

    confidence: float | None = None

    evidence: list[str] = field(default_factory=list)


# ---------------------------------------------------------------------
# Metadata
# ---------------------------------------------------------------------

@dataclass
class Metadata:
    """
    Provenance for this Variant.
    """

    uuid: UUID = field(default_factory=uuid4)

    created: datetime = field(default_factory=lambda: datetime.now(UTC))
    modified: datetime = field(default_factory=lambda: datetime.now(UTC))

    modelmaker_version: str = "0.3.0-dev"
    pipeline_version: str = "0.1"

    uniprot_release: str | None = None
    ensembl_release: str | None = None

    history: list[str] = field(default_factory=lambda: [
        "Variant created."
    ])


# ---------------------------------------------------------------------
# User Input
# ---------------------------------------------------------------------

@dataclass(frozen=True)
class VariantInput:
    """
    Original user-supplied biological hypothesis.

    This object is immutable.
    """

    species: str

    human_gene: str

    mutation: str

    transcript: str | None = None

    protein_accession: str | None = None

    notes: str | None = None


# ---------------------------------------------------------------------
# Placeholder knowledge containers
# ---------------------------------------------------------------------

@dataclass
class OrthologyInfo(EvidenceContainer):
    """
    Orthology relationship between the human gene and the model organism.
    """

    human_gene: str | None = None
    model_gene: str | None = None

    human_accession: str | None = None
    model_query: str | None = None

    identity: float | None = None
    similarity: float | None = None


@dataclass
class GeneFamily(EvidenceContainer):
    """
    Represents a gene family and its known members.

    Useful for paralogs, duplicated genes, and protein families
    where multiple potential model genes exist.
    """

    family_name: str | None = None

    human_gene: str | None = None

    members: list[str] = field(default_factory=list)

    notes: str | None = None

@dataclass
class ExpressionProfile(EvidenceContainer):
    """
    Tissue and developmental expression profile.
    """

    tissues: list[str] = field(default_factory=list)

    developmental_stages: list[str] = field(default_factory=list)

    source: str | None = None


@dataclass
class ModelCandidate(EvidenceContainer):
    """
    Best candidate experimental model for the human variant.
    """

    species: str | None = None

    gene: str | None = None

    protein: str | None = None

    transcript: str | None = None

    orthology_score: float | None = None

    residue_match: bool | None = None

    rationale: list[str] = field(default_factory=list)

    expression: ExpressionProfile | None = None


@dataclass
class ProteinInfo(EvidenceContainer):
    """
    Protein-level alignment and residue mapping.
    """

    aligned_human: str | None = None
    aligned_model: str | None = None

    alignment_score: float | None = None
    similarity: float | None = None
    gaps: float | None = None
    identity: float | None = None


@dataclass
class GenomeInfo(EvidenceContainer):
    """Genome-level mapping."""


@dataclass
class AnnotationInfo(EvidenceContainer):
    """Biological annotation."""


@dataclass
class EditingInfo(EvidenceContainer):
    """Genome editing recommendations."""


@dataclass
class ValidationInfo(EvidenceContainer):
    """Experimental validation strategy."""


@dataclass
class ReportInfo(EvidenceContainer):
    """Reporting and visualization."""


# ---------------------------------------------------------------------
# Variant
# ---------------------------------------------------------------------

@dataclass
class Variant:
    """
    Canonical biological representation of a disease variant.
    """

    input: VariantInput

    metadata: Metadata = field(default_factory=Metadata)

    orthology: OrthologyInfo | None = None
    gene_family: GeneFamily | None = None
    model: ModelCandidate | None = None
    protein: ProteinInfo | None = None
    genome: GenomeInfo | None = None
    annotation: AnnotationInfo | None = None
    editing: EditingInfo | None = None
    validation: ValidationInfo | None = None
    report: ReportInfo | None = None