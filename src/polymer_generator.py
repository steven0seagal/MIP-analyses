"""
Polymer generation utilities using PySoftK and RDKit
"""
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors
import pandas as pd
from typing import List, Optional
import tempfile
import os

# Import PySoftK modules
try:
    from pysoftk.linear_polymer.super_monomer import *
    from pysoftk.linear_polymer.linear_polymer import *
    from pysoftk.format_printers.format_mol import *
except ImportError:
    print("Warning: PySoftK not installed. Polymer generation methods will fail.")
    # Define placeholder functions if PySoftK is missing
    def make_linear_polymer(*args, **kwargs):
        raise NotImplementedError("PySoftK's linear polymer function is missing.")
    def make_alternating_copolymer(*args, **kwargs):
        raise NotImplementedError("PySoftK's alternating copolymer function is missing.")


class PolymerGenerator:
    """Generate polymers from monomers"""

    def __init__(self):
        self.polymers = []
        self.metadata = []

    def _add_terminal_groups(self, mol: Chem.Mol, terminal_group: str = "Br") -> Chem.Mol:
        """
        Add a single terminal group to the last carbon atom for polymerization.

        Args:
            mol: RDKit molecule
            terminal_group: Terminal group to add (default: "Br")

        Returns:
            Modified RDKit molecule with terminal group at the last C
        """
        mol_copy = Chem.RWMol(mol)

        # Find the last carbon atom (highest index carbon)
        carbon_indices = []
        for atom in mol_copy.GetAtoms():
            if atom.GetSymbol() == 'C':
                carbon_indices.append(atom.GetIdx())

        if not carbon_indices:
            raise ValueError("No carbon atoms found in molecule")

        # Get the last carbon
        last_c_idx = max(carbon_indices)
        last_c = mol_copy.GetAtomWithIdx(last_c_idx)

        # Check if it has available valence
        if last_c.GetTotalNumHs() > 0:
            # Add terminal group atom
            terminal_idx = mol_copy.AddAtom(Chem.Atom(terminal_group))
            mol_copy.AddBond(last_c_idx, terminal_idx, Chem.BondType.SINGLE)
        else:
            raise ValueError(f"Last carbon atom (idx {last_c_idx}) has no available hydrogen to replace")

        return mol_copy.GetMol()

    def _detect_terminal_groups(self, mol: Chem.Mol) -> str:
        """
        Detect terminal groups from molecule.

        Args:
            mol: RDKit molecule

        Returns:
            Symbol of detected terminal group
        """
        # Priority order for terminal group detection
        terminal_candidates = ['Br', 'I', 'Cl', 'F']

        for atom in mol.GetAtoms():
            # Check if atom is terminal (has only 1 bond)
            if atom.GetDegree() == 1:
                symbol = atom.GetSymbol()
                if symbol in terminal_candidates:
                    return symbol

        return None

    def generate_homopolymer(self, monomer_smiles: str, dp: int, name: Optional[str] = None,
                           terminal_group: str = "Br") -> Chem.Mol:
        """
        Generate linear homopolymer using PySoftK.

        Args:
            monomer_smiles: SMILES string of monomer (must contain terminal groups like Br)
            dp: Degree of polymerization
            name: Optional name for polymer
            terminal_group: Terminal group for polymerization (default: "Br")

        Returns:
            RDKit molecule object
        """
        # try:
            # Create RDKit molecule from SMILES
        mol = Chem.MolFromSmiles(monomer_smiles)

        if mol is None:
            raise ValueError(f"Invalid SMILES: {monomer_smiles}")

        # Check if monomer has terminal groups
        mol_str = Chem.MolToSmiles(mol)
        if terminal_group not in mol_str:
            raise ValueError(f"Monomer must contain terminal group '{terminal_group}'. "
                            f"Example: 'BrCC(Br)c1ccccc1' for styrene with Br terminals")

        # Create super-monomer (same monomer repeated)
        a = Sm(mol, mol, terminal_group)
        k = a.mon_to_poly()

        # Create linear polymer with desired degree of polymerization
        # Lp returns PySoftK Molecule object (not RDKit)
        pysoftk_polymer = Lp(k, terminal_group, dp, shift=1.0).linear_polymer()

        if pysoftk_polymer is None:
            raise ValueError(f"Failed to generate polymer from {monomer_smiles}")

        # Use PySoftK's Fmt to save as MOL file and read back with RDKit
        with tempfile.NamedTemporaryFile(mode='w', suffix='.mol', delete=False) as tmp:
            tmp_path = tmp.name

        # Use Fmt().mol_print() to save as MOL format (just writes file)
        Fmt(pysoftk_polymer).mol_print(tmp_path)
        print(tmp_path)
        # Read back with RDKit to get proper mol object
        polymer = Chem.MolFromMolFile(tmp_path, removeHs=False)
        os.unlink(tmp_path)

        if polymer is None:
            raise ValueError(f"Failed to convert PySoftK polymer to RDKit molecule")

        # Get SMILES for metadata
        polymer_smiles = Chem.MolToSmiles(polymer)

        # Store metadata
        metadata = {
            'name': name or f'Polymer_{len(self.polymers)}',
            'monomer': monomer_smiles,
            'dp': dp,
            'type': 'homopolymer',
            'smiles': polymer_smiles,
            'terminal_group': terminal_group
        }
        self.metadata.append(metadata)

        print(f"Generated {metadata['name']} (Homopolymer) with DP={dp}")
        return polymer

        # except Exception as e:
        #     print(f"Error generating homopolymer with PySoftK: {e}")
        #     import traceback
        #     traceback.print_exc()
        #     return None

    def generate_copolymer(self, monomer1: str, monomer2: str, type: str, dp: int,
                          name: Optional[str] = None, terminal_group: str = "Br") -> Chem.Mol:
        """
        Generate copolymer using PySoftK.

        Args:
            monomer1: First monomer SMILES (must contain terminal groups like Br)
            monomer2: Second monomer SMILES (must contain terminal groups like Br)
            type: Type of copolymer (e.g., 'alternating', 'block')
            dp: Total degree of polymerization (number of monomer units in the chain)
            name: Optional name
            terminal_group: Terminal group for polymerization (default: "Br")

        Returns:
            RDKit molecule object
        """
        try:
            # Create RDKit molecules from SMILES
            mol_1 = Chem.MolFromSmiles(monomer1)
            mol_2 = Chem.MolFromSmiles(monomer2)

            if mol_1 is None or mol_2 is None:
                raise ValueError(f"Invalid SMILES: {monomer1} or {monomer2}")

            # Check if monomers have terminal groups
            mol1_str = Chem.MolToSmiles(mol_1)
            mol2_str = Chem.MolToSmiles(mol_2)
            if terminal_group not in mol1_str or terminal_group not in mol2_str:
                raise ValueError(f"Monomers must contain terminal group '{terminal_group}'. "
                               f"Example: 'BrCC(Br)c1ccccc1' for styrene with Br terminals")

            # Create super-monomer from two different monomers
            a = Sm(mol_1, mol_2, terminal_group)
            k = a.mon_to_poly()

            # Create linear polymer with desired degree of polymerization
            # Lp returns PySoftK Molecule object (not RDKit)
            pysoftk_polymer = Lp(k, terminal_group, dp, shift=1.0).linear_polymer()

            if pysoftk_polymer is None:
                raise ValueError(f"Failed to generate copolymer from {monomer1} and {monomer2}")

            # Use PySoftK's Fmt to save as MOL file and read back with RDKit
            with tempfile.NamedTemporaryFile(mode='w', suffix='.mol', delete=False) as tmp:
                tmp_path = tmp.name

            # Use Fmt().mol_print() to save as MOL format (just writes file)
            Fmt(pysoftk_polymer).mol_print(tmp_path)

            # Read back with RDKit to get proper mol object
            polymer = Chem.MolFromMolFile(tmp_path, removeHs=False)
            os.unlink(tmp_path)

            if polymer is None:
                raise ValueError(f"Failed to convert PySoftK polymer to RDKit molecule")

            # Get SMILES for metadata
            polymer_smiles = Chem.MolToSmiles(polymer)

            # Store metadata
            metadata = {
                'name': name or f'Polymer_{len(self.polymers)}',
                'monomer1': monomer1,
                'monomer2': monomer2,
                'dp': dp,
                'type': f'{type} copolymer',
                'smiles': polymer_smiles,
                'terminal_group': terminal_group
            }
            self.metadata.append(metadata)

            print(f"Generated {metadata['name']} ({type} Copolymer) with DP={dp}")
            return polymer

        except Exception as e:
            print(f"Error generating copolymer with PySoftK: {e}")
            import traceback
            traceback.print_exc()
            return None

    def optimize_structure(self, mol: Chem.Mol, force_field: str = 'MMFF94') -> Chem.Mol:
        """
        Optimize polymer structure using force field
        ... (method remains the same)
        """
        try:
            mol_h = Chem.AddHs(mol)
            AllChem.EmbedMolecule(mol_h)

            if force_field == 'MMFF94':
                AllChem.MMFFOptimizeMolecule(mol_h)
            elif force_field == 'UFF':
                AllChem.UFFOptimizeMolecule(mol_h)

            return Chem.RemoveHs(mol_h)
        except Exception as e:
            print(f"Optimization error: {e}")
            return mol

    def save_library(self, filename: str):
        """Save polymer library to SDF file with embedded SMILES and metadata"""
        writer = Chem.SDWriter(filename)
        for i, mol in enumerate(self.polymers):
            if mol is not None:
                # Add metadata properties before writing
                meta = self.metadata[i]
                for key, value in meta.items():
                     mol.SetProp(key, str(value))

                # Ensure SMILES is stored in the SDF
                if 'smiles' in meta:
                    mol.SetProp('SMILES', meta['smiles'])

                writer.write(mol)
        writer.close()
        print(f"Saved {len(self.polymers)} polymers to {filename}")

    def generate_library(self, monomers: dict, dp_values: List[int]) -> pd.DataFrame:
        """
        Generate systematic polymer library
        ... (method adjusted to use self.metadata)
        """
        self.polymers = []
        self.metadata = []

        for name, smiles in monomers.items():
            for dp in dp_values:
                polymer = self.generate_homopolymer(smiles, dp, f"{name}_DP{dp}")
                if polymer:
                    self.polymers.append(polymer)
                    # Metadata is now stored in self.metadata within generate_homopolymer

        results = []
        for i, polymer in enumerate(self.polymers):
             meta = self.metadata[i]
             results.append({
                  'name': meta.get('name'),
                  'monomer': meta.get('monomer'),
                  'dp': meta.get('dp'),
                  'type': meta.get('type'),
                  'mw': Descriptors.MolWt(polymer)
             })

        return pd.DataFrame(results)