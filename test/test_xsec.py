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
        self.assertAlmostEqual(np.log10(self.csms.total_cross_section(1e9, 'nu', 'NC')),
            np.log10(1.7316076522013216e-30), 2)
        self.assertAlmostEqual(np.log10(self.csms.total_cross_section(1e15, 'nu', 'CC')),
            np.log10(1.7664654996566748e-24), 3)
        self.assertAlmostEqual(np.log10(self.dipole.total_cross_section(1e18, 'nu', 'NC')),
            np.log10(1.1564151708472636e-23), 3)
        self.assertAlmostEqual(np.log10(self.dipole.total_cross_section(1e15, 'nu', 'CC')),
            np.log10(1.908338488012135e-24), 3)

    def test_xs_raises(self):
        with self.assertRaises(AttributeError):
            self.dipole.total_cross_section(0.1, 'bar', 'NC')
        with self.assertRaises(AttributeError):
            self.csms.differential_cross_section(10., [1., 2., 3.], 'nubar', 'BC')

    def test_differential(self):
        #self.assertEqual(self.dipole.differential_cross_section(0.1, [1.0], 'nu', 'CC'),
            #0.)
        self.assertAlmostEqual(np.log10(self.dipole.differential_cross_section(1e15, np.asarray([0.1, 0.2, 0.3]), 'nu', 'CC')[1]),
            np.log10(2.776026967357417e-40), 3)
        self.assertAlmostEqual(np.log10(self.csms.differential_cross_section(1e15, np.asarray([0.1, 0.2, 0.3]), 'nu', 'CC')[1]),
            np.log10(6.7270713004890386e-40), 3)
        self.assertAlmostEqual(np.log10(self.csms.differential_cross_section(1e15, np.asarray([0.1, 0.2, 0.3]), 'nu', 'NC')[-1]),
            np.log10(2.7432728976770823e-40), 3)

if __name__ == '__main__':
    unittest.main()
