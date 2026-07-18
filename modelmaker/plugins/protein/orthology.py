"""
Protein orthology plugin.

Wraps the existing orthology/alignment engine and enriches a Variant.
"""
from modelmaker.align import align_proteins
from modelmaker.core.plugin import Plugin
from modelmaker.core.variant import (
    AlignmentInfo,
    GeneFamily,
    ModelCandidate,
    OrthologyInfo,
    ProteinInfo,
)

class OrthologyPlugin(Plugin):

    name = "Orthology"

    requires = []

    provides = ["orthology", "protein", "alignment", "gene_family", "model"]

    version = "0.1"

    def run(self, variant):
        print("Orthology plugin running")

        alignment = align_proteins(
            human_accession="P02708",
            zfish_query="gene:chrna1 AND organism_id:7955",
            human_gene="CHRNA1",
            zfish_gene="CHRNA1",
        )

        variant.orthology = OrthologyInfo(
            human_gene=alignment.human_gene,
            model_gene=alignment.zfish_gene,
            human_accession=alignment.human_accession,
            model_query=alignment.zfish_query,
            identity=alignment.identity,
            similarity=alignment.similarity,
            confidence=alignment.identity,
            evidence=[
                "Selected best UniProt ortholog by sequence identity."
            ],
        )

        variant.alignment = AlignmentInfo(
            human_protein_id=alignment.human_accession,
            model_protein_id=alignment.zfish_gene,
            human_sequence=alignment.human_sequence,
            model_sequence=alignment.zfish_sequence,
            aligned_human=alignment.aligned_human,
            aligned_model=alignment.aligned_zfish,
            identity=alignment.identity,
            similarity=alignment.similarity,
            alignment_score=alignment.score,
            gaps=alignment.gaps,
            confidence=alignment.identity,
            evidence=[
                "Protein alignment completed using legacy aligner."
            ],
        )

        variant.protein = ProteinInfo(
            protein_id=alignment.zfish_gene,
            gene_symbol=alignment.zfish_gene,
            species=variant.input.species,
            sequence=alignment.zfish_sequence,
            length=len(alignment.zfish_sequence),
            source="legacy aligner / UniProt",
            confidence=alignment.identity,
            evidence=[
                "Selected model protein sequence from legacy aligner."
            ],
        )

        variant.gene_family = GeneFamily(
            family_name=None,
            human_gene=alignment.human_gene,
            model_species=variant.input.species,
            members=[alignment.zfish_gene],
            preferred_member=alignment.zfish_gene,
            notes="Single ortholog identified.",
        )

        variant.model = ModelCandidate(
            species=variant.input.species,
            gene=alignment.zfish_gene,
            protein_id=alignment.zfish_gene,
            selection_score=alignment.identity,
            residue_conserved=True,
            expression_supported=None,
            confidence=alignment.identity,
            evidence=[
                "Best ortholog by sequence identity.",
            ],
            rationale=[
                "Best ortholog by sequence identity.",
            ],
        )

        variant.metadata.history.append(
            "Orthology, alignment, protein, gene family, and model candidate selected."
        )

        return variant