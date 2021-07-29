import unittest, os
from taurunner import Casino
from taurunner.main import *
from taurunner.body import Body
from taurunner.modules import make_initial_e
from taurunner.modules import make_initial_thetas

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

        for en in cls.atten_dict.keys():
            # set up MC
            TR_specs = {
                'energy': en,
                'theta': 0.,
                'flavor': 16, 
                'no_secondaries': True, 
                'nevents': num_sim,
                'seed': 5,
                'xs_model': 'dipole',
                'no_losses': True
                }
            rand = np.random.RandomState(TR_specs['seed'])
            TR_specs['rand'] = rand
            eini = make_initial_e(TR_specs, rand=rand)
            thetas = make_initial_thetas(TR_specs, rand=rand)
            body = Body(6.0, 500.0)
            tracks  = {theta:Chord(theta=theta, depth=0.) for theta in set(thetas)}
            xs = CrossSections(TR_specs['xs_model'])
            prop = make_propagator(body, xs_model=xs.model)
            sim = run_MC(eini, thetas, body, xs, tracks, TR_specs, prop)
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
        self.assertTrue(relative_diffs[self.sim_dict_ens[2]] < 0.2)


if __name__ == '__main__':
    unittest.main()
