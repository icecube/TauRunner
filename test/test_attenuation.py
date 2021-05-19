import unittest, os
from taurunner import Casino
from taurunner.main import run_MC
from taurunner.body import Body

import numpy as np

def find_nearest_ind(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx

class TestConservation(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        """ once before all tests """
        num_sim = 200
        ens = [1e5, 1e6, 1e7]
        test_path = os.path.dirname(os.path.realpath(__file__))
        cls.tabled_comparisons = np.load(f'{test_path}/.comparisons_for_tests/attenuation.npy')
        cls.tabled_ens = cls.tabled_comparisons[0]
        cls.tabled_attens = cls.tabled_comparisons[1]
        table_inds = [find_nearest_ind(cls.tabled_ens, en*1e9) for en in ens]
        cls.atten_dict = {cls.tabled_ens[ind]/1e9: cls.tabled_attens[ind] for ind in table_inds}
        cls.sim_dict = dict()
        body = Body(6.0, 500.0) # Construct earthsized body of constant density
        for en in cls.atten_dict.keys():
            sim = run_MC(num_sim, 1, theta=0.,
                                      debug=False, xs_model='dipole', 
                                      energy=en, body=body, 
                                      losses=False,
                                     )
            cls.sim_dict[en] = sim
        cls.num_sim = num_sim
        cls.sim_dict_ens = list(cls.sim_dict.keys())

    @classmethod
    def tearDownClass(cls):
        """ once after all tests """
        pass

    def setUp(self):
        """ before each test """
        pass

    def tearDown(self):
        """ after each test """
        pass

    def test_energies(self):
        en_key = self.sim_dict_ens[0]
        incoming_en = np.sum(self.sim_dict[en_key]['Eini'])
        outgoing_en = np.sum(self.sim_dict[en_key]['Eout'])
        self.assertTrue(incoming_en >= outgoing_en)

    def test_number_conservation(self):
        en_key = self.sim_dict_ens[-1]
        nutau_msk = self.sim_dict[en_key]['PDG_Encoding'] == 16
        outgoing_nutau = len(self.sim_dict[en_key]['Eini'][nutau_msk])
        self.assertTrue(outgoing_nutau <= self.num_sim)

    def test_attenuation(self):
        relative_diffs = {}
        for en in self.sim_dict.keys():
            msk = self.sim_dict[en]['nCC'] == 0
            msk *= self.sim_dict[en]['nNC'] == 0
            nutau_msk = self.sim_dict[en]['PDG_Encoding'] == 16
            noninteracting = len(self.sim_dict[en]['Eini'][msk*nutau_msk])
            survival_prob = float(noninteracting) / float(self.num_sim)
            diff = (survival_prob - self.atten_dict[en]) \
                / self.atten_dict[en]
            diff = np.abs(diff)
            relative_diffs[en] = diff
        self.assertTrue(relative_diffs[self.sim_dict_ens[0]] < 0.1)
        self.assertTrue(relative_diffs[self.sim_dict_ens[1]] < 0.1)
        self.assertTrue(relative_diffs[self.sim_dict_ens[2]] < 0.1)


if __name__ == '__main__':
    unittest.main()
