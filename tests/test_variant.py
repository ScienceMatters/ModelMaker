from modelmaker.core.variant import Variant
from modelmaker.core.variant import VariantInput


def test_variant_creation():

    variant = Variant(
        input=VariantInput(
            species="danio_rerio",
            human_gene="CHRNE",
            mutation="P121L",
        )
    )

    assert variant.input.human_gene == "CHRNE"

    assert variant.orthology is None

    assert variant.genome is None

    assert variant.metadata.uuid is not None

    assert len(variant.metadata.history) == 1