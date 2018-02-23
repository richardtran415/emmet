import logging
import os
from datetime import datetime

from monty.json import jsanitize
from monty.serialization import loadfn

from pymatgen.core.structure import Structure
from pymatgen.analysis.diffraction.xrd import XRDCalculator, WAVELENGTHS

from maggma.builder import Builder

__author__ = "Shyam Dwaraknath <shyamd@lbl.gov>"

module_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)))
default_xrd_settings = os.path.join(
    module_dir, "settings", "xrd.json")


class SurfaceBuilder(Builder):

    def __init__(self, materials, surface, query=None, **kwargs):
        """
        Calculates surface properties for materials

        Args:
            materials (Store): Store of materials documents
            surface (Store): Store of surface properties
            query (dict): dictionary to limit materials to be analyzed
        """

        self.materials = materials
        self.surface = surface
        self.query = query if query else {}

        super().__init__(sources=[materials],
                         targets=[surface],
                         **kwargs)

    def get_items(self):
        """
        Gets all materials that need a new XRD

        Returns:
            generator of materials to calculate xrd
        """

        self.logger.info("Surface Builder Started")

        self.logger.info("Setting indexes")
        self.ensure_indexes()

        # All relevant materials that have been updated since diffraction props
        # were last calculated
        q = dict(self.query)
        q.update(self.materials.lu_filter(self.surface))
        mats = list(self.materials.distinct(self.materials.key, q))
        self.logger.info(
            "Found {} new materials for surface data".format(len(mats)))
        for m in mats:
            yield self.materials.query(properties=[self.materials.key, "structure"],
                                       criteria={self.materials.key: m}).limit(1)[0]

    def process_item(self, item):
        """
        Calculates surface properties for the structures

        Args:
            item (dict): a dict with a material_id and a structure

        Returns:
            dict: a surface dict
        """
        self.logger.debug("Calculating surface properties for {}".format(
            item[self.materials.key]))

        struct = Structure.from_dict(item['structure'])
        surf_doc

        return surf_doc

    def update_targets(self, items):
        """
        Inserts the new task_types into the task_types collection

        Args:
            items ([[dict]]): a list of list of surface dictionaries to update
        """
        items = list(filter(None, items))

        if len(items) > 0:
            self.logger.info("Updating {} surface docs".format(len(items)))
            self.surface.update(docs=items)
        else:
            self.logger.info("No items to update")
