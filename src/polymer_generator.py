"""
Polymer generation utilities using PySoftK and RDKit
"""
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors
import pandas as pd
from typing import List, Optional

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
        Add terminal groups to a molecule for polymerization.
        Finds terminal carbon atoms and adds halogen groups.

        Args:
            mol: RDKit molecule
            terminal_group: Terminal group to add (default: "Br")

        Returns:
            Modified RDKit molecule with terminal groups
        """
        # Find terminal carbons (degree 1, only H neighbors)
        terminal_carbons = []
        for atom in mol.GetAtoms():
            if atom.GetSymbol() == 'C' and atom.GetDegree() <= 2:
                # Check if it has available valence
                h_count = atom.GetTotalNumHs()
                if h_count > 0:
                    terminal_carbons.append(atom.GetIdx())

        if len(terminal_carbons) < 2:
            raise ValueError(f"Molecule needs at least 2 positions to add terminal groups")

        # Add terminal groups to first and last terminal carbons
        # We'll use SMARTS replacement
        mol_copy = Chem.RWMol(mol)

        # Sort to get first and last
        terminal_carbons = sorted(terminal_carbons)[:2]

        # Add Br to these positions by replacing H
        for idx in reversed(terminal_carbons):  # reversed to maintain indices
            atom = mol_copy.GetAtomWithIdx(idx)
            if atom.GetTotalNumHs() > 0:
                # Add Br atom
                br_idx = mol_copy.AddAtom(Chem.Atom(terminal_group))
                mol_copy.AddBond(idx, br_idx, Chem.BondType.SINGLE)
                # Remove one implicit H
                atom.SetNumExplicitHs(max(0, atom.GetNumExplicitHs()))

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
            monomer_smiles: SMILES string of monomer
            dp: Degree of polymerization
            name: Optional name for polymer
            terminal_group: Terminal group for polymerization (default: "Br")

        Returns:
            RDKit molecule object
        """
        try:
            # Create RDKit molecule from SMILES
            mol = Chem.MolFromSmiles(monomer_smiles)

            if mol is None:
                raise ValueError(f"Invalid SMILES: {monomer_smiles}")

            # Detect if terminal groups already exist
            detected_group = self._detect_terminal_groups(mol)

            if detected_group is None:
                # Add terminal groups to monomer
                print(f"Adding terminal groups ({terminal_group}) to monomer")
                mol = self._add_terminal_groups(mol, terminal_group)
                print(f"Modified monomer SMILES: {Chem.MolToSmiles(mol)}")
            else:
                terminal_group = detected_group
                print(f"Using existing terminal group: {terminal_group}")

            # Create super-monomer (same monomer repeated)
            a = Sm(mol, mol, terminal_group)
            k = a.mon_to_poly()

            # Create linear polymer with desired degree of polymerization
            polymer = Lp(k, terminal_group, dp, shift=1.0).linear_polymer()

            if polymer is None:
                raise ValueError(f"Failed to generate polymer from {monomer_smiles}")

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

        except Exception as e:
            print(f"Error generating homopolymer with PySoftK: {e}")
            import traceback
            traceback.print_exc()
            return None

    def generate_copolymer(self, monomer1: str, monomer2: str, type: str, dp: int,
                          name: Optional[str] = None, terminal_group: str = "Br") -> Chem.Mol:
        """
        Generate copolymer using PySoftK.

        Args:
            monomer1: First monomer SMILES
            monomer2: Second monomer SMILES
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

            # Detect or add terminal groups
            detected_group = self._detect_terminal_groups(mol_1)

            if detected_group is None:
                print(f"Adding terminal groups ({terminal_group}) to monomers")
                mol_1 = self._add_terminal_groups(mol_1, terminal_group)
                mol_2 = self._add_terminal_groups(mol_2, terminal_group)
                print(f"Modified monomer1 SMILES: {Chem.MolToSmiles(mol_1)}")
                print(f"Modified monomer2 SMILES: {Chem.MolToSmiles(mol_2)}")
            else:
                terminal_group = detected_group
                print(f"Using existing terminal group: {terminal_group}")

            # Create super-monomer from two different monomers
            a = Sm(mol_1, mol_2, terminal_group)
            k = a.mon_to_poly()

            # Create linear polymer with desired degree of polymerization
            polymer = Lp(k, terminal_group, dp, shift=1.0).linear_polymer()

            if polymer is None:
                raise ValueError(f"Failed to generate copolymer from {monomer1} and {monomer2}")

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