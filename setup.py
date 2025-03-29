from setuptools import setup, find_packages

setup(
    name="magnumopus",
    version="0.1.0",
    packages=find_packages(),
    install_requires=[
        "numpy>=1.21.0",
        "pandas>=1.3.0",
        "biopython>=1.79",
        "setuptools>=45.0.0",
        "wheel>=0.37.0",
    ],
    author="Sanya Badole",
    author_email="sanyabadole1@gmail.com",
    description="A DNA sequence analysis tool for PCR primer analysis and sequence alignment",
    long_description=open("README.md").read(),
    long_description_content_type="text/markdown",
    url="https://github.com/sanyabadole/DNA-Sequence-Analysis",
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Intended Audience :: Science/Research",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
    python_requires=">=3.8",
    entry_points={
        "console_scripts": [
            "magnumopus=magnumopus.magop:main",
        ],
    },
) 
