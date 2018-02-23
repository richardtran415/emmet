import unittest
import os

from emmet.vasp.builders.elastic import *
from maggma.stores import MongoStore
from maggma.runner import Runner

from monty.serialization import loadfn

__author__ = "Richard Tran"
__version__ = "0.1"
__maintainer__ = "Richard Tran"
__email__ = "rit001@eng.ucsd.edu"
__date__ = "02/23/18"

# Test with provided tasks for max Miller index
# of 1 for Li (bcc), Si (diamond) and Co (hcp).
module_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)))
test_tasks = os.path.join(module_dir, "..","..","..", "..", "test_files",
                          "vasp", "surface_tasks.json")

class SurfaceBuilderTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        # Set up test db, set up mpsft, etc.
        cls.test_tasks = MongoStore("test_emmet", "tasks")
        cls.test_tasks.connect()
        docs = loadfn(test_tasks, cls=None)
        cls.test_tasks.update(docs)
        cls.test_surface = MongoStore("test_emmet", "surface")

        # Generate test materials collection
        cls.test_materials = MongoStore("test_emmet", "materials")
        cls.test_materials.connect()
        # opt_docs = cls.test_tasks.query(["output.structure", "formula_pretty"],
        #                                 {"task_label": "structure optimization"})
        # mat_docs = [{"material_id": "mp-{}".format(n),
        #              "structure": opt_doc['output']['structure'],
        #              "pretty_formula": opt_doc['formula_pretty']}
        #             for n, opt_doc in enumerate(opt_docs)]
        # cls.test_materials.update(mat_docs, key='material_id', update_lu=False)
