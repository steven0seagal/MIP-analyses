#!/bin/bash

echo "Installing external polymer tools..."

# Create directory for external tools
mkdir -p external_tools
cd external_tools

# Install PySoftK
echo ""
echo "1/2: Installing PySoftK..."
if [ ! -d "PySoftK" ]; then
    git clone https://github.com/SoftMatterTheoryGroup/PySoftK.git
    cd PySoftK
    pip install -e .
    cd ..
    echo "✓ PySoftK installed"
else
    echo "✓ PySoftK already exists"
fi

# Install m2p
echo ""
echo "2/2: Installing m2p..."
if [ ! -d "m2p" ]; then
    git clone https://github.com/NREL/m2p.git
    cd m2p
    pip install -e .
    cd ..
    echo "✓ m2p installed"
else
    echo "✓ m2p already exists"
fi

cd ..

echo ""
echo "External tools installation complete!"
