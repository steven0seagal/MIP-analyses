#!/bin/bash

echo "Installing MIP Polymer Toolkit..."

if command -v conda &> /dev/null; then
    conda env create -f environment.yml
    echo "Environment created. Activate with: conda activate mip-toolkit"
else
    python3 -m venv venv
    source venv/bin/activate
    pip install -r requirements.txt
    echo "Environment created. Activate with: source venv/bin/activate"
    cd external_tools && git clone https://github.com/alejandrosantanabonilla/pysoftk.git && pip install .
fi

pip install -e .
echo "Setup complete!"
