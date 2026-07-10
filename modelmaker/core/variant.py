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

    Every inference should include:
      - confidence
      - evidence supporting the inference
      - provenance (filled by future plugins)
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
    """Orthology relationships."""


@dataclass
class ProteinInfo(EvidenceContainer):
    """Protein-level mapping."""


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
    protein: ProteinInfo | None = None
    genome: GenomeInfo | None = None
    annotation: AnnotationInfo | None = None
    editing: EditingInfo | None = None
    validation: ValidationInfo | None = None
    report: ReportInfo | None = None