"""
Central biological object used throughout ModelMaker.

The Variant contains biological state only.

It never performs computations.

Every plugin enriches the Variant with additional biological knowledge.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from datetime import datetime
from typing import Any


# ---------------------------------------------------------------------
# Metadata
# ---------------------------------------------------------------------

@dataclass
class Metadata:
    modelmaker_version: str = "0.3.0-dev"
    created: datetime = field(default_factory=datetime.utcnow)
    modified: datetime = field(default_factory=datetime.utcnow)

    uniprot_release: str | None = None
    ensembl_release: str | None = None

    notes: list[str] = field(default_factory=list)


# ---------------------------------------------------------------------
# User Input
# ---------------------------------------------------------------------

@dataclass
class VariantInput:

    species: str

    gene: str

    mutation: str

    transcript: str | None = None

    protein_accession: str | None = None

    notes: str | None = None


# ---------------------------------------------------------------------
# Generic data containers
# ---------------------------------------------------------------------

@dataclass
class OrthologyInfo:

    confidence: float | None = None

    evidence: list[str] = field(default_factory=list)

    data: dict[str, Any] = field(default_factory=dict)


@dataclass
class ProteinInfo:

    confidence: float | None = None

    evidence: list[str] = field(default_factory=list)

    data: dict[str, Any] = field(default_factory=dict)


@dataclass
class GenomeInfo:

    confidence: float | None = None

    evidence: list[str] = field(default_factory=list)

    data: dict[str, Any] = field(default_factory=dict)


@dataclass
class AnnotationInfo:

    confidence: float | None = None

    evidence: list[str] = field(default_factory=list)

    data: dict[str, Any] = field(default_factory=dict)


@dataclass
class EditingInfo:

    confidence: float | None = None

    evidence: list[str] = field(default_factory=list)

    data: dict[str, Any] = field(default_factory=dict)


@dataclass
class ValidationInfo:

    confidence: float | None = None

    evidence: list[str] = field(default_factory=list)

    data: dict[str, Any] = field(default_factory=dict)


@dataclass
class ReportInfo:

    confidence: float | None = None

    evidence: list[str] = field(default_factory=list)

    data: dict[str, Any] = field(default_factory=dict)


# ---------------------------------------------------------------------
# Variant
# ---------------------------------------------------------------------

@dataclass
class Variant:

    input: VariantInput

    metadata: Metadata = field(default_factory=Metadata)

    orthology: OrthologyInfo = field(default_factory=OrthologyInfo)

    protein: ProteinInfo = field(default_factory=ProteinInfo)

    genome: GenomeInfo = field(default_factory=GenomeInfo)

    annotation: AnnotationInfo = field(default_factory=AnnotationInfo)

    editing: EditingInfo = field(default_factory=EditingInfo)

    validation: ValidationInfo = field(default_factory=ValidationInfo)

    report: ReportInfo = field(default_factory=ReportInfo)