import unittest, os
from taurunner import Casino
from taurunner.main import *
from taurunner.utils import make_initial_e
from taurunner.utils import make_initial_thetas
from taurunner.body import Body
from taurunner.utils.make_tracks import make_tracks

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
        rand = np.random.RandomState(2)
        eini = make_initial_e(num_sim, en)
        xs = CrossSections('dipole')
        thetas = make_initial_thetas(num_sim, 0.)
        body = Body(6.0, 500.)
        tracks = make_tracks(thetas)
        sim = run_MC(
            eini, 
            thetas, 
            body, 
            xs, 
            tracks, 
            None,
            no_secondaries=False,
            no_losses=True
        )
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
