import unittest

import numpy as np
from taurunner.modules import construct_body
import taurunner
import taurunner.body as body

from taurunner.modules import units

class TestBodyMethods(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        """ once before all tests """
        rad_km = 500.0
        dens = 6.
        tst_bod = construct_body({'body': 6., 'radius': rad_km})
        tst_earth = construct_body({'body': 'earth', 'water': 0.})
        tst_sun = construct_body({'body': 'sun'})
        cls.bod_param = {'rad_km': rad_km, 'dens': dens}
        cls.body = tst_bod
        cls.earth = tst_earth
        cls.sun = tst_sun

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

    def test_body_equalities(self):
        self.assertEqual(self.body.radius, self.bod_param['rad_km']*units.km)
        self.assertEqual(len(self.body.layer_boundaries), 2)
        self.assertEqual(self.body.layer_boundaries[0], 0.)
        self.assertEqual(self.body.layer_boundaries[-1], 1.)
        self.assertEqual(self.body.get_density(0.5),
            self.bod_param['dens']*(units.gr / units.cm**3))

    def test_body_errors(self):
        with self.assertRaises(IndexError):
            self.body.get_density(1.5)
        with self.assertRaises(TypeError):
            self.body.radius()

    def test_earth(self):
        self.assertEqual(self.earth.radius/units.km, 
            6368.0)
        dens_arr = np.array([self.earth.get_density(rr) for rr in np.r_[0.:1.0:20j]]) / \
            (units.gr / units.cm**3)
        ref_arr = np.array([13.0885, 13.06401773, 12.99057091, 12.86815956, 12.10241034,
            11.89591806, 11.64506629, 11.34501924, 10.99094114, 10.57799621,
            10.10134867,  5.46234005,  5.29539984,  5.12674874,  4.95369186,
            4.77353433,  4.58358127,  4.38113778,  3.50463684,  2.6       ])
        self.assertTrue((np.around(dens_arr, 3) == np.around(ref_arr, 3)).all())
        self.assertTrue( (np.around(self.earth.layer_boundaries * self.earth.radius / units.km, 3) == \
            np.around(np.asarray([0., 1221., 3480., 5701., 5771., 5971., 6151., 6346.6, 6356., 6368.]), 3)).all())

    def test_HZ_sun(self):
        self.assertAlmostEqual(self.sun.get_density(0.5), 5.78148105641597e+18,
            4)
        self.assertAlmostEqual(self.sun.get_edensity(0.5), 1.0287332840597812e-14,
            4)
        self.assertAlmostEqual(self.sun.radius / units.km, 696300.0,
            3)
        self.assertAlmostEqual(self.sun.get_edensity(1), 1.5130182132307767e-21,
            5)

    def test_bad_body(self):
        with self.assertRaises(RuntimeError):
            tmp = body.Body([1., 2., 3.], 3.)
        with self.assertRaises(RuntimeError):
            tmp = body.Body([1., 2., 3.], 3., 
                layer_boundaries=[0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6])
    
    def test_total_mass(self):
        import logging
        logging.getLogger('scipy').setLevel(logging.ERROR)
        mass = body.check_total_mass(self.earth)[0]
        self.assertAlmostEqual(mass / units.kg, 5.970273732299387e+24, 5)

if __name__ == '__main__':
    unittest.main()