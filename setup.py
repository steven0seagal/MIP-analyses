from setuptools import setup, find_packages

setup(
    name="mip-polymer-toolkit",
    version="0.1.0",
    author="Your Name",
    author_email="your.email@example.com",
    description="Computational toolkit for MIP polymer analysis",
    packages=find_packages(),
    python_requires=">=3.8",
    install_requires=[
        "rdkit>=2023.3.1",
        "numpy>=1.24.0",
        "pandas>=2.0.0",
        "matplotlib>=3.7.0",
        "seaborn>=0.12.0",
        "scikit-learn>=1.3.0",
    ],
)
