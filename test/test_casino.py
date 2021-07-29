import unittest, os
from taurunner import Casino
from taurunner.particle import Particle
from taurunner.cross_sections import CrossSections
import taurunner

import numpy as np

class TestCasinoHelpers(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        """ once before all tests """
        dip_xs = CrossSections('dipole')
        csms_xs = CrossSections('CSMS')
        dipole_particle = Particle(15, 1e12, 0., 0., 1, dip_xs, None, None, True, False)
        csms_particle = Particle(15, 1e12, 0., 0., 1, csms_xs, None, None, False, False)
        nutau = Particle(16, 1e12, 0., 0., 1, dip_xs, None, None, True, False)
        mu = Particle(13, 1e12, 0., 0., 1, dip_xs, None, None, False, False) 
        cls.dipole = dipole_particle
        cls.csms = csms_particle
        cls.nutau = nutau
        cls.mu = mu
        

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

    def test_particle_getters(self):
        self.assertEqual(self.dipole.GetParticleId(), 15)
        self.assertEqual(self.csms.GetLifetime(), self.dipole.GetLifetime())
        self.assertAlmostEqual(self.csms.GetLifetime(), 441.67, 3)
        self.assertAlmostEqual(self.csms.GetMass(), 1776000000.0, 4)
        self.assertAlmostEqual(self.nutau.GetInteractionDepth('CC'), 
            5.545643041570203e+34, 5)
        self.assertAlmostEqual(self.nutau.GetInteractionProbability(1e33, 'CC'),
            0.017870567135038762, 5)

    def test_particle_decays(self):
        self.mu.Decay()
        self.assertTrue(~self.mu.survived)

    def test_particle_raises(self):
        with self.assertRaises(ValueError):
            self.dipole.GetInteractionDepth('CC')
        with self.assertRaises(ValueError):
            self.nutau.Decay()
        with self.assertRaises(ValueError):
            self.dipole.Interact('CC')


if __name__ == '__main__':
    unittest.main()
