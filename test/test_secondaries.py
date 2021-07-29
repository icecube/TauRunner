import unittest, os
from taurunner import Casino
from taurunner.main import *
from taurunner.modules import construct_body
from taurunner.modules import make_initial_e
from taurunner.modules import make_initial_thetas

import numpy as np

def find_nearest_ind(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx

class TestSecondaries(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        """ once before all tests """
        num_sim = 200
        en = 1e6
        TR_specs = {
            'flavor': 16,
            'no_secondaries': False,
            'nevents': num_sim,
            'seed': 2,
            'energy': en,
            'theta': 0.,
            'xs_model': 'dipole',
            'no_losses': True
            }
        rand = np.random.RandomState(TR_specs['seed'])
        TR_specs['rand'] = rand
        body = construct_body({'body': 6.0, 'radius': 500.})
        xs = CrossSections('dipole')
        eini = make_initial_e(TR_specs, rand=rand)
        thetas = make_initial_thetas(TR_specs)
        tracks  = {theta:Chord(theta=theta, depth=0.) for theta in set(thetas)}
        sim = run_MC(
            eini, 
            thetas, 
            body, 
            xs, 
            tracks, 
            TR_specs,
            None)
        cls.sim = sim
        cls.num_sim = num_sim
        cls.en = en

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
        incoming_en = self.en * self.num_sim * 1e9
        outgoing_en = np.sum(self.sim['Eout'])
        self.assertTrue(incoming_en >= outgoing_en)

    def test_broken_number_conservation(self):
        ''' With secondaries, there should be more outgoing particles
        than there are incoming particles'''
        self.assertTrue(len(self.sim['Eout']) >= self.num_sim)

    def test_number_of_secondaries(self):
        '''Rough numbers ro make sure the number of secondaries
        compares well to expectations'''
        nutaus = self.sim['PDG_Encoding'] == 16
        secondaries = ~nutaus
        self.assertTrue(len(self.sim[nutaus]) <= self.num_sim)
        self.assertTrue(len(self.sim[nutaus]) > 0.9*self.num_sim)
        self.assertTrue(len(self.sim[~nutaus][self.sim[~nutaus]['Eout'] != 0.]) < 0.20 * self.num_sim)

if __name__ == '__main__':
    unittest.main()
