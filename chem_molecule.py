from molecule import Molecule
import time
from rdkit.Chem import AllChem as Chem


class ChemMolecule(Molecule):

    """
    Molecules are conceptually just combinations of Atoms, described by some character string symbol.
    The symbol has no particular meaning except to an IChemistry. For example, an IChemistry based on RDKit would expect symbols in SMILES form.

    Any transformations of Molecules is the responsibility of a IChemistry.
    """

    def __init__(self, symbol):
        """
        Initialise a new ChemMolecule represented by character string in SMILES notation.

        :param symbol: SMILES string
        """

        if not isinstance(symbol, basestring):
            raise TypeError

        mol = Chem.MolFromSmiles(symbol)
        if mol is None:
            raise ValueError("SMILES:{} is either not legal or not sensible".format(symbol))

        mol = Chem.AddHs(mol)   # MolFromSmiles doesn't add Hs even if SMILES contains explicit Hs

        self.mass = sum([atom.GetMass() for atom in mol.GetAtoms()])
        self._symbol = Chem.MolToSmiles(mol)  # will include Hs
        self._id = "{}.{}".format(id(self), time.clock())

    def get_symbol(self):
        return self._symbol

    def __str__(self):
        return self._symbol
