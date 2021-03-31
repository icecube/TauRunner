import unittest
from taurunner import Casino
from taurunner.main import propagate_neutrinos

import numpy as np

# include number conservation
# include attenuation check
# include energy conservation check 
class TestConservation(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        """ once before all tests """
        num_sim = 1000
        ens = [1e5, 1e6, 1e7]
        cls.sim_dict = dict()
        for en in ens:
            sim = propagate_neutrinos(num_sim, 1, theta=60., 
                debug=False, xs_model='dipole', energy=en)
            cls.sim_dict[en] = sim
        cls.num_sim = num_sim
        cls.atten_dict = {1e5: 0.63, 1e6: 0.22, 1e7: 0.012}

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
        incoming_en = 0.
        outgoing_en = 0.
        for ptype in self.sim_dict[1e6]:
            incoming_en += np.sum(ptype['Eini'])
            outgoing_en += np.sum(ptype['Eout'])
        self.assertTrue(incoming_en >= outgoing_en)

    def test_number_conservation(self):
        outgoing_n = 0
        for ptype in self.sim_dict[1e6]:
            outgoing_n += len(ptype['Eini'])
        self.assertTrue(outgoing_n == self.num_sim)

    def test_attenuation(self):
        relative_diffs = {}
        for en in self.sim_dict.keys():
            noninteracting = 0
            for ptype in self.sim_dict[en]:
                msk = ptype['nCC'] == 0
                msk *= ptype['nNC'] == 0
                noninteracting += len(ptype['Eini'][msk])
            survival_prob = float(noninteracting) / float(self.num_sim)
            diff = (survival_prob - self.atten_dict[en]) \
                / self.atten_dict[en]
            diff = np.abs(diff)
            relative_diffs[en] = diff
            print(en, survival_prob, diff)
        self.assertTrue(relative_diffs[1e5] < 0.1)
        self.assertTrue(relative_diffs[1e6] < 0.1)
        # increase threshold for 1e7 because low statistics
        self.assertTrue(relative_diffs[1e7] < 0.5)


if __name__ == '__main__':
    unittest.main()