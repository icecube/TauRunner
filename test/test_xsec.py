import unittest
import numpy as np

from taurunner.cross_sections import CrossSections

class TestXSecMethods(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        """ once before all tests """
        dipole = CrossSections('dipole')
        csms = CrossSections('CSMS')
        cls.dipole = dipole
        cls.csms = csms

    def test_total(self):
        self.assertAlmostEqual(np.log10(self.csms.TotalNeutrinoCrossSection(1e9)),
            np.log10(1.7316076522013216e-30), 2)
        self.assertAlmostEqual(np.log10(self.csms.TotalNeutrinoCrossSection(1e15, interaction='CC')),
            np.log10(1.7664654996566748e-24), 3)
        self.assertAlmostEqual(np.log10(self.dipole.TotalNeutrinoCrossSection(1e18)),
            np.log10(1.1564151708472636e-23), 3)
        self.assertAlmostEqual(np.log10(self.dipole.TotalNeutrinoCrossSection(1e15, interaction='CC')),
            np.log10(1.908338488012135e-24), 3)

    def test_xs_raises(self):
        with self.assertRaises(ValueError):
            self.dipole.TotalNeutrinoCrossSection(0.1)
        with self.assertRaises(ValueError):
            self.csms.DifferentialOutGoingLeptonDistribution(10., [1., 2., 3.], 'BC')

    def test_differential(self):
        self.assertEqual(self.dipole.DifferentialOutGoingLeptonDistribution(0.1, [1.0], 'CC')[0],
            0.)
        self.assertAlmostEqual(np.log10(self.csms.DifferentialOutGoingLeptonDistribution(30, [1., 2., 3.], 'CC')[1]),
            np.log10(6.965933260169092e-30), 3)
        self.assertAlmostEqual(np.log10(self.csms.DifferentialOutGoingLeptonDistribution(30., [1., 2., 3.], 'NC')[-1]),
            np.log10(2.244095228466293e-30), 3)

if __name__ == '__main__':
    unittest.main()