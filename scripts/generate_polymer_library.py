import sys
sys.path.insert(0, 'src')

from polymer_generator import PolymerGenerator
import pandas as pd

# Define monomers for library
monomers = {
    'styrene': 'C=Cc1ccccc1',
    'ethylene': 'C=C',
    'methyl_acrylate': 'C=CC(=O)OC',
    'vinyl_acetate': 'C=COC(=O)C',
    'propylene': 'CC=C',
    'vinyl_chloride': 'C=CCl'
}

# Degree of polymerization values
dp_values = [10, 20, 50, 100]

# Initialize generator
gen = PolymerGenerator()

print(f"\nGenerating polymers for {len(monomers)} monomers with {len(dp_values)} DP values...")
print(f"Total polymers to generate: {len(monomers) * len(dp_values)}")
print("-" * 60)

# Generate library
df = gen.generate_library(monomers, dp_values)

# Save results
output_file = 'data/polymers/library.sdf'
gen.save_library(output_file)

# Save metadata
df.to_csv('data/polymers/library_metadata.csv', index=False)

print("\n" + "=" * 60)
print("Library generation complete!")
print(f"  Polymers: {len(gen.polymers)}")
print(f"  Output: {output_file}")
print(f"  Metadata: data/polymers/library_metadata.csv")
print("=" * 60)