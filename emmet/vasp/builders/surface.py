import logging
import os
from datetime import datetime

from monty.json import jsanitize
from monty.serialization import loadfn

from pymatgen.core.structure import Structure
from pymatgen import MPRester
from pymatgen.analysis.wulff import WulffShape

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
    """
    Format of each entry in surface_properties collection:

        {"material_id": str,
         "polymorph": int,
         "weighted_surface_energy_EV_PER_ANG2": float,
         "weigthed_work_function": float,
         "e_above_hull": float,
         "pretty_formula": str,
         "weighted_surface_energy": float,
         "anisotropy": float,
         "shape_factor": float,
         "spacegroup": {"symbol": str,
                        "number": int},
         "surfaces":
           [
               {"miller_index": list,
                "tasks": {"OUC": int, "slab": int},
                "surface_energy_EV_PER_ANG2": float,
                "surface_energy": float,
                "work_function": float,
                "is_reconstructed": bool,
                "area_fraction": float,
                "structure": str(cif)
               },
               {
               },
               {
               },
               etc ...
           ]
        }

    Format of each entry in wulff collection:

        {"material_id": str,
         "polymorph": int,
         "pretty_formula": str,
         "spacegroup": {"symbol": str,
                        "number": int},
         "thumbnail": binary,
         "hi_res_images":
           [
               {"miller_index": list,
                "image": binary
               },
               {
               },
               {
               },
               etc ...
           ]
        }
    """

    def __init__(self, materials, surface, query=None, MAPI_KEY=None, **kwargs):
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
        self.mprester = MPRester(api_key=MAPI_KEY)

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
        f = self.materials.lu_filter(surfprops)
        q.update(f)

        mats = list(self.materials.distinct(materials.key, q))
        self.logger.info(
            "Found {} new slab calculations for {} materials".format(len(mats)))

        for m in mats:
            slab_task = self.materials.query(criteria={self.materials.key: m}).limit(1)[0]
            ouc_criteria = {"structure_type": "oriented_unit_cell",
                            "miller_index": slab_task["miller_index"],
                            "reconstruction": slab_task["reconstruction"]}
            ouc_task = self.materials.query(criteria=ouc_criteria).limit(1)[0]

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
        surface_entry["initial_structure"] = slab_task["initial_structure"]
        surface_entry["work_function"] = l - slab_outputs["efermi"]

        # Get the surface properties for this material
        el = ouc[0].species_string
        surfprops = self.get_material_surfprops(slab_task["material_id"], el)

        # Determine whether or not to include this termination in the collection
        surfprops = self.check_termination_stability(surfprops, surface_entry)

        # Get wulff shape related qunatities
        surfprops = self.get_wulff_surfprops(surfprops)

        return surfprops

    def get_material_surfprops(self, mpid, el):
        """
        Gets existing surface properties of a material from the collection.
            If non-existent make an empty entry (no slab surface properties)
            before inserting surface properties.
        """

        surf_props = self.surface.find_one({"material_id": mpid})
        if not surf_props:
            surf_props = {}
            entries = self.mprester.get_entries(el, property_data=["e_above_hull",
                                                                   "material_id",
                                                                   "spacegroup"])
            sort_entries = sorted(entries, key=lambda entry: entry.data["e_above_hull"])
            for i, entry in enumerate(sort_entries):
                if entry.data["material_id"] == mpid:
                    surf_props["polymorph"] = i
                    surf_props["e_above_hull"] = entry.data["e_above_hull"]
                    surf_props["material_id"] =  mpid
                    surf_props["pretty_formula"] = el
                    surf_props["spacegroup"] = {"symbol": entry.data["spacegroup"]["symbol"],
                                                "number": entry.data["spacegroup"]["number"]}
                    surf_props["surfaces"] = []

                    self.logger.info("No surface properties for {}-{}. Creating new doc.".format(el, mpid))
                    self.surface.insert(surf_props)

                break

        return surf_props

    def check_termination_stability(self, surfprops, propdoc):
        """
        For a particular surface of Miller index hkl, check to see if the
            current termination in question is the most stable (relative to
            whatever entries reconstructed or not exists in the collection).

        Three possibilities:
            1) there are no facet entries as of yet
            2) only one entry exists (its unreconstructed)
            3) two entries exists (its reconstructed)

        :return:
        """

        surfaces = []
        facetprops = surfprops["surfaces"]
        se = propdoc["surface_energy"]
        # Do entries for this specific facet exist already?
        if not facetprops:
            # If not, then just append the current facet entry
            surfaces.append(propdoc)
        elif len(facetprops) > 1:
            # That means one of these entries must be reconstructed
            resurf = [surf for surf in facetprops if surf['is_reconstructed']][0]
            surf = [surf for surf in facetprops if not surf['is_reconstructed']][0]
            if not propdoc["is_reconstructed"] and \
                    all([se < surf["surface_energy"],
                         se < resurf["surface_energy"]]):
                # then replace both entries, reconstruction won't happen and
                # the current unreconstructed termination is too unstable
                surfaces.append(propdoc)
            elif not propdoc["is_reconstructed"] \
                    and propdoc["surface_energy"] < surf["surface_energy"]:
                # It is only more stable than the unreconstructed
                # entry, replace that entry only
                surfaces.append(propdoc)
                surfaces.append(resurf)
            elif propdoc["is_reconstructed"] and \
                            propdoc["surface_energy"] < resurf["surface_energy"]:
                # it is a more stable reconstruction, append it along with the unreconstructed
                surfaces.append(propdoc)
                surfaces.append(surf)
            else:
                # do nothing
                surfaces.extend(facetprops)
        else:
            # That means no reconstruction as of yet
            if propdoc["surface_energy"] < facetprops[0]["surface_energy"]:
                # replace the current doc
                surfaces.append(propdoc)
            else:
                # do not do anything
                surfaces.extend(facetprops)

        surfprops["surfaces"] = surfaces

        return surfprops

    def get_wulff_shape(self, surfprops):
        """
        Get the Wulff shape
        :param surfprops:
        :return:
        """

        e_surf_list = [surf["surface_energy"] for surf in surfprops["surfaces"]]
        miller_list = [surf["miller_index"] for surf in surfprops["surfaces"]]
        conv = self.materials.find_one({"material_id": surfprops["material_id"],
                                        "structure_type": "conventional_unit_cell"})
        ucell = Structure.from_str(conv["initial_structure"], "cif")
        return WulffShape(ucell.lattice, miller_list, e_surf_list)

    def get_wulff_surfprops(self, surfprops):
        """
        Get surface properties determined from the Wulff shape
        :return:
        """

        wulffshape = self.get_wulff_shape(surfprops)
        se = wulffshape.weighted_surface_energy

        area_frac_dict = wulffshape.area_fraction_dict
        weighted_work_function = 0
        for surface in surfprops["surfaces"]:
            hkl = tuple(surface["miller_index"])
            f = area_frac_dict[hkl]
            surface["area_fraction"] = f
            weighted_work_function += f*surface["work_function"]

        surfprops["weighted_work_function"] = weighted_work_function / wulffshape.surface_area
        surfprops["shape_factor"] = wulffshape.shape_factor
        surfprops["weighted_surface_energy"] = se
        surfprops["weighted_surface_energy_EV_PER_ANG2"] = se/EV_PER_ANG2_TO_JOULES_PER_M2
        surfprops["surface_anisotropy"] = wulffshape.anisotropy

        return surfprops

    def update_targets(self, surfprops):
        """
        Make sure to insert the slab_cell task_id as an
        individual item to make querying for new entries easier?
        Also make sure this particular slab is the most stable
        one before inserting.
        """

        self.logger.info("Updating {}-{} surface properties".format(surfprops["pretty_formula"],
                                                                    surfprops["material_id"]))
        self.surface.update_one({"material_id": surfprops["material_id"]},
                                {"$set": surfprops})

