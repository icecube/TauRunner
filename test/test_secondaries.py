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
        en = 1e6
        body = Body(6.0, 500.0) 
        sim = run_MC(num_sim, 2, theta=0.,
            debug=False, xs_model='dipole', 
            energy=en, body=body, 
            losses=False,
            with_secondaries=True
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
        # self.assertTrue(len(self.sim[~nutaus][self.sim[~nutaus]['Eout'] != 0.]) < 0.20 * self.num_sim)


if __name__ == '__main__':
    unittest.main()