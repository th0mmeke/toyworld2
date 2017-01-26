import unittest
import cycle_utilities


class TestCycleUtilities(unittest.TestCase):
    def test_flatten(self):
        self.assertListEqual([1, 2, 3, 4], cycle_utilities.flatten([[1, 2], [3, 4]]))

    def test_get_molecules(self):
        self.assertItemsEqual(['1a', '1b'], cycle_utilities.get_molecules('1a+1b>'))
        self.assertItemsEqual(['1a', '1b'], cycle_utilities.get_molecules('>1a+1b'))
        self.assertItemsEqual(['a'], cycle_utilities.get_molecules('a>'))
        self.assertItemsEqual(['1a', '1b', '2c'], cycle_utilities.get_molecules('1a+2c+1b>'))

    def test_get_products(self):
        self.assertEqual({'2c'}, cycle_utilities.get_products(['1a+1b>', '>2c']))
        self.assertEqual({'2c', '4a'}, cycle_utilities.get_products(['1a+1b>', '>2c+4a']))
        self.assertEqual({'2c', '1d', '4a'}, cycle_utilities.get_products(['1a+1b>', '>2c+4a', '>1d']))

    def test_get_reactants(self):
        self.assertEqual({'2c'}, cycle_utilities.get_reactants(['>1a+1b', '2c>']))
        self.assertEqual({'2c', '4a'}, cycle_utilities.get_reactants(['>1a+1b', '2c+4a>']))
        self.assertEqual({'2c', '1d', '4a'}, cycle_utilities.get_reactants(['>1a+1b', '2c+4a>', '1d>']))

    def test_get_molecules_in_cycle(self):
        self.assertSetEqual({'1a', '1b', '2c'}, cycle_utilities.get_molecules_in_cycle(['>1a+1b', '2c>']))
        self.assertSetEqual({'2c', '1d', '4a', '1a', '1b'},
                            cycle_utilities.get_molecules_in_cycle(['>1a+1b', '2c+4a>', '1d>']))

    def test_map_id_to_smiles(self):
        smiles = { '1': 'a',
                   '2': 'b',
                   '3': 'c',
                   '4': 'd',
                   '5': 'e',
                   }
        molecular_cycles = ['1', '1+2>', '>3+4+5', '3']
        species_cycle = ['a', 'a+b>', '>c+d+e', 'c']
        self.assertListEqual(species_cycle, cycle_utilities.map_id_to_smiles(molecular_cycles, smiles))

    def test_discover_species(self):
        smiles = { '1': 'a',
                   '2': 'b',
                   '3': 'c',
                   '4': 'd',
                   '5': 'e',
                   }
        molecular_cycles = [{'cycle': ['1', '1+2>', '>3+4+5', '3']}, {'cycle': ['1', '1+2>', '>3+4+5', '3']}]
        expected = {frozenset(['a', 'c', 'a+b>', '>c+d+e']): [['1', '1+2>', '>3+4+5', '3'], ['1', '1+2>', '>3+4+5', '3']]}
        self.assertEqual(expected, cycle_utilities.discover_species(molecular_cycles, smiles, length=0))

    def test_identify_clusters(self):
        molecular_cycles = [['1', '1+2>', '>3+4+5', '3'], ['3+6>', '>7']]
        actual = cycle_utilities.identify_clusters(molecular_cycles)
        expected = [[['1', '1+2>', '>3+4+5', '3'], ['3+6>', '>7']]]
        self.assertListEqual(actual, expected)

        molecular_cycles = [['1', '1+2>', '>3+4+5', '3'], ['8+6>', '>7']]
        actual = cycle_utilities.identify_clusters(molecular_cycles)
        expected = []
        self.assertListEqual(actual, expected)

        molecular_cycles = [['1', '1+2>', '>3+4+5', '3'], ['4>', '>9']]
        actual = cycle_utilities.identify_clusters(molecular_cycles)
        expected = [[['1', '1+2>', '>3+4+5', '3'], ['4>', '>9']]]
        self.assertListEqual(actual, expected)

    def test_discover_autocatalysis(self):
        clusters = [[['1', '1+2>', '>3+4+5', '3'], ['3+6>', '>7']]]
        actual = cycle_utilities.discover_autocatalysis(clusters)
        expected = []
        self.assertListEqual(actual, expected)

        clusters = [[['1', '1+2>', '>3+4+5', '3'], ['4>', '>9']]]
        actual = cycle_utilities.discover_autocatalysis(clusters)
        expected = []
        self.assertListEqual(actual, expected)

        clusters = [[['1', '1+2>', '>3+4+5', '3'], ['3+6>', '>7'], ['4+9>', '>10']]]
        actual = cycle_utilities.discover_autocatalysis(clusters)
        expected = [[[], ['3'], ['4']]]  # linkage molecules in downstream cycles
        self.assertListEqual(actual, expected)

        clusters = [[['1', '1+2>', '>3+4+5', '3'], ['3+6>', '>7'], ['4+9>', '>10']], [['11', '11+12>', '>13+14+15'], ['13+16>', '>17'], ['14+19>', '>20']]]
        actual = cycle_utilities.discover_autocatalysis(clusters)
        expected = [[[], ['3'], ['4']], [[], ['13'], ['14']]]
        self.assertListEqual(actual, expected)

        clusters = [[['1', '1+2>', '>3+4+5', '3'], ['3+6>', '>7'], ['4+9>', '>10']], [['13+16>', '>17'], ['11', '11+12>', '>13+14+15'], ['14+19>', '>20']]]
        actual = cycle_utilities.discover_autocatalysis(clusters)
        print(actual)
        expected = [[[], ['3'], ['4']], [['13'], [], ['14']]]
        self.assertListEqual(actual, expected)

        clusters = [[['1', '1+2>', '>3+4+5'], ['3+10>', '>7'], ['4+3>', '>10']]]
        actual = cycle_utilities.discover_autocatalysis(clusters)
        expected = [[[], ['3'], ['3', '4']]]  # linkage molecules in downstream cycles
        self.assertListEqual(actual, expected)
