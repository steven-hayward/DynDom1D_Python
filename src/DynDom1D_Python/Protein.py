import gemmi
import numpy as np
"""
Gemmi follows a hierarchy:
Structure -> Model -> Chain -> Residue -> Atom
"""


class Protein:
    def __init__(self, input_path, protein_name, chain, atoms_to_use):
        self.input_path = input_path
        self.name = protein_name
        self.file_path = f"{input_path}/{self.name}.pdb"
        self.chain_param: str = chain  # The chain specified in parameter input for the program
        self.atoms_to_use = atoms_to_use
        # The index of the residues in the ResidueSpan
        self.utilised_residues_indices = []
        self.slide_window_residues_indices = None

    def get_structure(self):
        return gemmi.read_structure(self.file_path, format=gemmi.CoorFormat.Pdb)

    def get_model(self):
        structure = gemmi.read_structure(self.file_path, format=gemmi.CoorFormat.Pdb)
        return structure[0]

    def get_chain(self):
        structure = gemmi.read_structure(self.file_path, format=gemmi.CoorFormat.Pdb)
        return structure[0][self.chain_param]

    def get_polymer(self):
        structure = gemmi.read_structure(self.file_path, format=gemmi.CoorFormat.Pdb)
        return structure[0][self.chain_param].get_polymer()

    def get_utilised_residues(self):
        chain: gemmi.ResidueSpan = self.get_polymer()
        slide_window_chain = gemmi.Chain(self.chain_param)
        for i in range(len(self.utilised_residues_indices)):
            index = self.utilised_residues_indices[i]
            slide_window_chain.add_residue(chain[index])
        return slide_window_chain.get_polymer()

    def get_slide_window_residues(self):
        """
        Get the residues from only the sliding window result.
        :return:
        """
        chain: gemmi.ResidueSpan = self.get_polymer()
        slide_window_chain = gemmi.Chain(self.chain_param)
        for i in range(self.slide_window_residues_indices[0], self.slide_window_residues_indices[1]):
            index = self.utilised_residues_indices[i]
            slide_window_chain.add_residue(chain[index])
        return slide_window_chain.get_polymer()

    def print_chain(self):
        atoms = []
        res_span = self.get_polymer()
        for i in self.utilised_residues_indices:
            has_N = False
            has_CA = False
            has_C = False
            res_atoms = []

            for atom in res_span[i]:
                if atom.name == "N" and not has_N:
                    res_atoms.append(atom)
                elif atom.name == "CA" and not has_CA:
                    res_atoms.append(atom)
                elif atom.name == "C" and not has_C:
                    res_atoms.append(atom)
                if len(res_atoms) >= 3:
                    break

        backbone_atoms = np.asarray(atoms)
        print(backbone_atoms)
        print(f"{self.name}({self.chain_param}) - {backbone_atoms.shape}")
        print(f"{self.name}({self.chain_param}) - {backbone_atoms}")

