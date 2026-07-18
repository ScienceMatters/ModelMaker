from modelmaker.core.plugin import Plugin
from modelmaker.core.variant import AnnotationInfo


class ProteinAnnotationPlugin(Plugin):

    name = "ProteinAnnotation"

    requires = ["protein"]

    provides = ["annotation"]

    version = "0.1"

    def run(self, variant):

        print("Protein annotation plugin running")

        variant.annotation = AnnotationInfo(
            protein_name=variant.protein.gene_symbol,
            source="placeholder"
        )

        variant.metadata.history.append(
            "Protein annotation completed."
        )

        return variant