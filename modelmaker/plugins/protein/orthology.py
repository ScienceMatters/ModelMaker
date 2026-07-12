"""
Protein orthology plugin.

Wraps the existing orthology/alignment engine and enriches a Variant.
"""
from modelmaker.core.variant import OrthologyInfo
from modelmaker.align import align_proteins
from modelmaker.core.plugin import Plugin


class OrthologyPlugin(Plugin):

    name = "Orthology"

    requires = []

    provides = ["orthology", "protein"]

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

        variant.metadata.history.append(
            "Orthology completed."
        )

        return variant