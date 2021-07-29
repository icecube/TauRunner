import unittest, os, copy
import numpy as np

import taurunner
from taurunner.main import run_MC
from taurunner.track import Chord
from taurunner.cross_sections import CrossSections
from taurunner.modules import construct_body
from taurunner.modules import make_initial_e
from taurunner.modules import make_initial_thetas
from taurunner.modules import make_propagator

class TestMainMethods(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        """ once before all tests """
        path = os.path.dirname(os.path.realpath(__file__))
        num_sim = 30
        # set up MC
        TR_specs_base = {
            'flavor': 16,
            'no_secondaries': True,
            'nevents': num_sim,
            'seed': 5,
            'xs_model': 'dipole',
            'no_losses': True
            }
        rand = np.random.RandomState(TR_specs_base['seed'])
        TR_specs_base['rand'] = rand
        body = construct_body({'body': 'earth', 'water': 0.})
        xs = CrossSections('dipole')
        prop = make_propagator(body, xs_model=xs.model)

        # GZK test
        TR_specs = copy.deepcopy(TR_specs_base)
        TR_specs['energy'] = path + '/../taurunner/gzk_cdf_phi_spline.npy'
        TR_specs['theta'] = 'range'
        TR_specs['th_max'] = 90.
        TR_specs['th_min'] = 0.
        gzk_eini = make_initial_e(TR_specs, rand=rand)
        gzk_thetas = make_initial_thetas(TR_specs)
        gzk_tracks  = {theta:Chord(theta=theta, depth=0.) for theta in set(gzk_thetas)}
        gzk = run_MC(
            gzk_eini, 
            gzk_thetas, 
            body, 
            xs, 
            gzk_tracks, 
            TR_specs, 
            prop)
        cls.gzk = gzk

        # spec_1
        TR_specs = copy.deepcopy(TR_specs_base)
        TR_specs['energy'] = -1.
        TR_specs['e_max'] = 1e8
        TR_specs['e_min'] = 1e5
        TR_specs['theta'] = 'range'
        TR_specs['th_max'] = 90.
        TR_specs['th_min'] = 0.
        spec_1_eini = make_initial_e(TR_specs, rand=rand)
        spec_1_thetas = make_initial_thetas(TR_specs)
        spec_1_tracks  = {theta:Chord(theta=theta, depth=0.) for theta in set(spec_1_thetas)}
        spec_1 = run_MC(
            spec_1_eini, 
            spec_1_thetas, 
            body, 
            xs, 
            spec_1_tracks, 
            TR_specs, 
            prop)
        cls.spec_1 = spec_1

        # spec 2
        TR_specs = copy.deepcopy(TR_specs_base)
        TR_specs['energy'] = -2.
        TR_specs['e_max'] = 1e8
        TR_specs['e_min'] = 1e5
        TR_specs['theta'] = 30.
        TR_specs['no_secondaries'] = False
        spec_2_eini = make_initial_e(TR_specs, rand=rand)
        spec_2_thetas = make_initial_thetas(TR_specs)
        spec_2_tracks  = {theta:Chord(theta=theta, depth=0.) for theta in set(spec_2_thetas)}
        spec_2 = run_MC(
            spec_2_eini, 
            spec_2_thetas, 
            body, 
            xs, 
            spec_2_tracks, 
            TR_specs,
            prop)
        cls.spec_2 = spec_2

        # mono
        TR_specs = copy.deepcopy(TR_specs_base)
        TR_specs['energy'] = 1e8
        TR_specs['theta'] = 'range'
        TR_specs['th_max'] = 90.
        TR_specs['th_min'] = 0.
        mono_eini = make_initial_e(TR_specs, rand=rand)
        mono_thetas = make_initial_thetas(TR_specs)
        mono_tracks  = {theta:Chord(theta=theta, depth=0.) for theta in set(mono_thetas)}
        mono = run_MC(
            mono_eini, 
            mono_thetas, 
            body, 
            xs, 
            mono_tracks, 
            TR_specs,
            prop)
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
