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
            np.log10(1.6697662677233126e-30), 2)
        self.assertAlmostEqual(np.log10(self.csms.TotalNeutrinoCrossSection(1e6, interaction='CC')),
            np.log10(1.9004075425072155e-36), 3)
        self.assertAlmostEqual(np.log10(self.dipole.TotalNeutrinoCrossSection(1e6)),
            np.log10(2.4515161246528807e-34), 3)
        self.assertAlmostEqual(np.log10(self.dipole.TotalNeutrinoCrossSection(1e6, interaction='CC')),
            np.log10(6.966542578303724e-34), 3)

    def test_xs_raises(self):
        with self.assertRaises(ValueError):
            self.dipole.TotalNeutrinoCrossSection(0.1)
        with self.assertRaises(ValueError):
            self.csms.DifferentialOutGoingLeptonDistribution(10., [1., 2., 3.], 'BC')

    def test_differential(self):
        self.assertEqual(self.dipole.DifferentialOutGoingLeptonDistribution(0.1, [1.0], 'CC')[0],
            0.)
        self.assertAlmostEqual(np.log10(self.csms.DifferentialOutGoingLeptonDistribution(30., [1., 2., 3.], 'CC')[1]),
            np.log10(6.34642047e-39), 3)
        self.assertAlmostEqual(np.log10(self.csms.DifferentialOutGoingLeptonDistribution(30., [1., 2., 3.], 'NC')[-1]),
            np.log10(2.02971510e-39), 3)

if __name__ == '__main__':
    unittest.main()