import unittest, os, copy
import numpy as np

import taurunner
from taurunner.main import *
from taurunner.cross_sections import CrossSections
from taurunner.body import Body
from taurunner.utils import make_initial_e
from taurunner.utils import make_initial_thetas
from taurunner.utils import make_propagator
from taurunner.utils.make_tracks import make_tracks

class TestMainMethods(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        """ once before all tests """
        path = os.path.dirname(os.path.realpath(__file__))

        num_sim = 30
        rand = np.random.RandomState(5)
        body = construct_earth(layers=[])
        xs = CrossSections('dipole')
        prop = make_propagator(16, body, xs_model=xs.model)

        # GZK test
        gzk_thetas = make_initial_thetas(num_sim, (0., 90.), rand=rand)
        gzk_eini = make_initial_e(
            num_sim, 
            path + '/../taurunner/resources/gzk_cdf_phi_spline.npy', 
            rand=rand
        )
        gzk = run_MC(
            gzk_eini, 
            gzk_thetas, 
            body, 
            xs, 
            prop,
            no_losses=True,
            no_secondaries=True
        )
        cls.gzk = gzk

        # spec_1
        spec_1_eini = make_initial_e(num_sim, -1., e_min=1e5, e_max=1e8, rand=rand)
        spec_1_thetas = make_initial_thetas(num_sim, (0., 90.), rand=rand)
        spec_1 = run_MC(
            spec_1_eini, 
            spec_1_thetas, 
            body, 
            xs, 
            prop,
            no_losses=True,
            no_secondaries=True
        )
        cls.spec_1 = spec_1

        # spec 2
        spec_2_eini = make_initial_e(num_sim, -2., e_min=1e5, e_max=1e8, rand=rand)
        spec_2_thetas = make_initial_thetas(num_sim, 30., rand=rand)
        spec_2 = run_MC(
            spec_2_eini, 
            spec_2_thetas, 
            body, 
            xs,
            prop,
            no_losses=True,
            no_secondaries=False
        )
        cls.spec_2 = spec_2

        # mono
        mono_eini = make_initial_e(num_sim, 1e8, rand=rand)
        mono_thetas = make_initial_thetas(num_sim, (0., 90.), rand=rand)
        mono = run_MC(
            mono_eini, 
            mono_thetas, 
            body, 
            xs,
            prop,
            no_losses=True,
            no_secondaries=True
        )
        cls.mono = mono

    def test_main_things(self):
        # print(np.median(self.gzk['Eini']))
        pass

    def test_main_raises(self):
        pass

    def test_initial_gzk_energies(self):
        pass

    def test_initial_spectral_sampling(self):
        pass

    def test_outgoing_gzk(self):
        pass

    # TODO: Reimplement all of the test_main_raises, but with new formatting,
    # Check some features from a gzk spectrum, including median simulated
    # and outgoing disrtibution. Check the expected features from a powerlaw
    # sampling

    # def test_main_raises(self):
    #     with self.assertRaises(RuntimeError):
    #         a = run_MC(None, 1)
    #     with self.assertRaises(RuntimeError):
    #         b = run_MC(100, 2, gzk=None, theta=None, energy=1e6)
    #     with self.assertRaises(RuntimeError):
    #         b = run_MC(100, 2, theta=30., energy=None, spectrum=None)
    #     with self.assertRaises(ValueError):
    #         b = run_MC(100, 2, spectrum=-2., theta=100., e_range="1e5 1e8")
    #     with self.assertRaises(ValueError):
    #         b = run_MC(100, 2, energy=1e7, theta=100., spectrum=None)
        # with self.assertRaises(ValueError):
        #     self.dipole.TotalNeutrinoCrossSection(0.1)
        # with self.assertRaises(ValueError):
        #     self.csms.DifferentialOutGoingLeptonDistribution(10., [1., 2., 3.], 'BC')

if __name__ == '__main__':
    unittest.main()
