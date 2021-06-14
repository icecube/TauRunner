import unittest
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
        self.assertAlmostEqual(self.csms.TotalNeutrinoCrossSection(1e9),
            1.1272239923709325e-23, 2)
        self.assertAlmostEqual(self.csms.TotalNeutrinoCrossSection(1e6, interaction='CC'),
            1.7661531211119893e-24, 3)
        self.assertAlmostEqual(self.dipole.TotalNeutrinoCrossSection(1e6),
            2.4515161246528807e-34, 3)
        self.assertAlmostEqual(self.dipole.TotalNeutrinoCrossSection(1e6, interaction='CC'),
            6.966542578303724e-34, 3)

if __name__ == '__main__':
    unittest.main()