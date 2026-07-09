"""
ModelMaker pipeline runner.
"""

from __future__ import annotations


class Pipeline:

    def __init__(self, plugins):

        self.plugins = plugins


    def run(self, variant):

        for plugin in self.plugins:

            variant = plugin.run(variant)

        return variant