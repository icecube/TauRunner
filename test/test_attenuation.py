import unittest
# from taurunner import Casino

# include number conservation
# include attenuation check
# include energy conservation check
class TestConservation(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        """ THIS IS RUN ONCE BEFORE ANY TESTS """
        pass

    @classmethod
    def tearDownClass(cls):
        """ once after all tests """
        pass

    def setUp(self):
        """THIS IS RUN BEFORE EACH INDIVIDUAL TEST"""
        pass

    def tearDown(self):
        """This is run after each test"""
        pass

    def test_upper(self):
        """Any class method that starts with test_ will be run and
        counts as a new test"""
        self.assertEqual('foo'.upper(), 'FOO')

    def test_isupper(self):
        self.assertTrue('FOO'.isupper())
        self.assertFalse('Foo'.isupper())

    def test_split(self):
        s = 'hello world'
        self.assertEqual(s.split(), ['hello', 'world'])
        # check that s.split fails when the separator is not a string
        with self.assertRaises(TypeError):
            s.split(2)

if __name__ == '__main__':
    unittest.main()
