import logging
import os
from datetime import datetime

from monty.json import jsanitize
from monty.serialization import loadfn

from pymatgen.core.structure import Structure

from maggma.builder import Builder

__author__ = "Richard Tran"
__version__ = "0.1"
__maintainer__ = "Richard Tran"
__email__ = "rit001@eng.ucsd.edu"
__date__ = "02/23/18"

module_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)))
default_xrd_settings = os.path.join(
    module_dir, "settings", "surface.json")


class SurfaceBuilder(Builder):
    def __init__(self, materials, surface, query=None, **kwargs):
        """
        Calculates surface properties for materials

        Args:
            materials (Store): Store of task documents for the raw
                calculations (via parsing with atomate's drone)
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
        Get any new slab_cell entry and their corresponding oriented_unit_cell
            tasks to be processed into the surface store yet
        """

        self.logger.info("Surface Builder Started")

        # I have no idea what this does
        # self.logger.info("Setting indexes")
        # self.ensure_indexes()

        q = {dict(self.query)}  # query criteria
        q["state"] = "successful"  # only find successfully completed calculations
        # only look for slab_cell calculations, because we need
        # a slab cell to actually process surface properties
        q["structure_type"] = "slab_cell"
        f = materials.lu_filter(surfprops)
        q.update(f)

        mats = list(materials.distinct(materials.key, q))
        self.logger.info(
            "Found {} new slab calculations for {} materials".format(len(mats)))

        for m in mats:
            slab_task = materials.query(criteria={materials.key: m}).limit(1)[0]
            ouc_criteria = {"structure_type": "oriented_unit_cell",
                            "miller_index": slab_task["miller_index"],
                            "reconstruction": slab_task["reconstruction"]}
            ouc_task = materials.query(criteria=ouc_criteria).limit(1)[0]

            yield slab_task, ouc_task

    def process_item(self, slab_task, ouc_task):
        surface_entry = {}
        hkl = slab_task["miller_index"]
        slab_outputs = slab_task["calcs_reversed"][0]["output"]
        ouc_outputs = ouc_task["calcs_reversed"][0]["output"]
        s = Structure.from_dict(outputs["structure"])
        slab_entry = SlabEntry(s, outputs["energy"], hkl)
        ouc = Structure.from_dict(ouc_outputs["structure"])
        ouc_entry = ComputedStructureEntry(ouc, ouc_outputs["energy"])
        se = slab_entry.surface_energy(ouc_entry)
        l = max(slab_outputs["locpot"]["2"])

        surface_entry["surface_energy"] = se * EV_PER_ANG2_TO_JOULES_PER_M2
        surface_entry["surface_energy_EV_PER_ANG2"] = se
        surface_entry["miller_index"] = hkl
        surface_entry["tasks"] = {'OUC': ouc_task["task_id"],
                                  'slab': slab_task["task_id"]}
        surface_entry["is_reconstructed"] = slab_task["reconstruction"]
        surface_entry["structure"] = s.to("cif")
        surface_entry["ouc_structure"] = ouc.to("cif")
        surface_entry["initial_structure"] = \
            Structure.from_dict(slab_task["initial_structure"])
        surface_entry["work_function"] = l - slab_outputs["efermi"]

        return surface_entry

    def update_targets(self):
        """
        Make sure to insert the slab_cell task_id as an
        individual item to make querying for new entries easier?
        Also make sure this particular slab is the most stable
        one before inserting.
        """
        pass
