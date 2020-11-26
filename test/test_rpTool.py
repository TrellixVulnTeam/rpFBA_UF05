import unittest
import os
import sys

sys.path.insert(0, '..')

import rpTool as rpFBA
#WARNING: Need to copy a version of rpSBML locally
import rpSBML

class TestRPTool(unittest.TestCase):

    """
    @classmethod
    def setUpClass(self):
    """

    def test_runFBA(self):
        rpsbml = rpSBML.rpSBML('test', path=os.path.join('data', 'merged.xml'))
        rpfba = rpFBA.rpFBA(rpsbml)
        obj_value, status = rpfba.runFBA('RP1_sink')
        self.assertAlmostEqual(obj_value, 9.230769230769237)
        self.assertTrue(status)
        #make sure that the results are written to the file
        all_json = rpsbml.genJSON()
        self.assertAlmostEqual(all_json['pathway']['brsynth']['fba_obj_RP1_sink']['value'], 9.230769230769237)

    def test_runFractionReaction(self):
        rpsbml = rpSBML.rpSBML('test', path=os.path.join('data', 'merged.xml'))
        rpfba = rpFBA.rpFBA(rpsbml)
        obj_value, status = rpfba.runFractionReaction('biomass', 1.0, 'RP1_sink', 1.0)
        self.assertAlmostEqual(obj_value, 2.3076923076923888)
        self.assertTrue(status)
        #make sure that the results are written to the file
        all_json = rpsbml.genJSON()
        self.assertAlmostEqual(all_json['pathway']['brsynth']['fba_obj_RP1_sink__restricted_biomass']['value'], 2.3076923076923888)
        self.assertAlmostEqual(all_json['pathway']['brsynth']['fba_obj_biomass']['value'], 3.6794124272706443)

    def test_runParsimoniousFBA(self):
        rpsbml = rpSBML.rpSBML('test', path=os.path.join('data', 'merged.xml'))
        rpfba = rpFBA.rpFBA(rpsbml)
        obj_value, status = rpfba.runParsimoniousFBA('RP1_sink')
        self.assertAlmostEqual(obj_value, 859.3846153846168)
        self.assertTrue(status)
        #make sure that the results are written to the file
        all_json = rpsbml.genJSON()
        self.assertAlmostEqual(all_json['pathway']['brsynth']['fba_obj_RP1_sink']['value'], 859.3846153846168)
