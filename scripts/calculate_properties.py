import sys
sys.path.insert(0, 'src')

from descriptor_calculator import PolymerDescriptorCalculator
from utils import load_polymers
import pandas as pd

# Load polymer library
print("Loading polymers...")
try:
    mols = load_polymers('data/polymers/library.sdf')
    print(f"✓ Loaded {len(mols)} polymers")
except FileNotFoundError:
    print("❌ Polymer library not found!")
    print("Run: bash scripts/generate_polymer_library.sh")
    sys.exit(1)

# Calculate descriptors
print("\nCalculating descriptors...")
calc = PolymerDescriptorCalculator()
df = calc.calculate_batch(mols)

# Save results
output_file = 'data/properties/descriptors.csv'
df.to_csv(output_file, index=False)

print("\n" + "=" * 60)
print("Property calculation complete!")
print(f"  Polymers analyzed: {len(df)}")
print(f"  Descriptors calculated: {len(df.columns) - 2}")  # Minus Name and SMILES
print(f"  Output: {output_file}")
print("=" * 60)

# Display summary statistics
print("\nSummary Statistics:")
print(df.describe().round(2))