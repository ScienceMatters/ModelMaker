"""
Base class for every ModelMaker module.
"""

from abc import ABC
from abc import abstractmethod


class Plugin(ABC):

    name = "Plugin"

    requires = []

    provides = []

    version = "0.1"


    @abstractmethod
    def run(self, variant):

        """Enrich a Variant."""

        pass