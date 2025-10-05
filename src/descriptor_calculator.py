"""
Comprehensive polymer descriptor calculator
"""
from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors, Lipinski
import pandas as pd
import numpy as np
from typing import List, Dict


class PolymerDescriptorCalculator:
    """Calculate molecular descriptors for polymers"""

    def __init__(self):
        self.descriptor_functions = {
            'MW': Descriptors.MolWt,
            'LogP': Descriptors.MolLogP,
            'TPSA': Descriptors.TPSA,
            'NumRotatableBonds': Descriptors.NumRotatableBonds,
            'NumHDonors': Descriptors.NumHDonors,
            'NumHAcceptors': Descriptors.NumHAcceptors,
            'FractionCsp3': Descriptors.FractionCSP3,
            'NumAromaticRings': Descriptors.NumAromaticRings,
            'NumSaturatedRings': Descriptors.NumSaturatedRings,
            'NumAliphaticRings': Descriptors.NumAliphaticRings,
        }

    def calculate_all(self, mol: Chem.Mol) -> Dict[str, float]:
        """
        Calculate all descriptors for a molecule

        Args:
            mol: RDKit molecule

        Returns:
            Dictionary of descriptor values
        """
        results = {}
        for name, func in self.descriptor_functions.items():
            try:
                results[name] = func(mol)
            except Exception as e:
                print(f"Error calculating {name}: {e}")
                results[name] = None
        return results

    def calculate_batch(self, mols: List[Chem.Mol], names: List[str] = None) -> pd.DataFrame:
        """
        Calculate descriptors for multiple molecules

        Args:
            mols: List of RDKit molecules
            names: Optional list of molecule names

        Returns:
            DataFrame with all descriptors
        """
        data = []
        for i, mol in enumerate(mols):
            if mol is None:
                continue
            row = self.calculate_all(mol)
            row['Name'] = names[i] if names else f"Polymer_{i}"
            row['SMILES'] = Chem.MolToSmiles(mol)
            data.append(row)

        df = pd.DataFrame(data)
        # Reorder columns
        cols = ['Name', 'SMILES'] + [c for c in df.columns if c not in ['Name', 'SMILES']]
        return df[cols]

    def calculate_mw_distribution(self, polymers: List[Chem.Mol]) -> Dict[str, float]:
        """
        Calculate molecular weight distribution statistics

        Args:
            polymers: List of polymer molecules

        Returns:
            Dictionary of MW statistics
        """
        mws = [Descriptors.MolWt(p) for p in polymers if p is not None]
        return {
            'mean': np.mean(mws),
            'std': np.std(mws),
            'min': np.min(mws),
            'max': np.max(mws),
            'median': np.median(mws)
        }

    def add_custom_descriptor(self, name: str, func):
        """Add custom descriptor function"""
        self.descriptor_functions[name] = func
        print(f"Added custom descriptor: {name}")
