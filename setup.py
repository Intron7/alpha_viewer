
# Always prefer setuptools over distutils
from setuptools import setup, find_packages
import pathlib

here = pathlib.Path(__file__).parent.resolve()

# Get the long description from the README file
long_description = (here / "README.md").read_text(encoding="utf-8")

# Arguments marked as "Required" below must be included for upload to PyPI.
# Fields marked as "Optional" may be commented out.

with open('requirements.txt') as f:
    required = f.read().splitlines()

setup(
    name="alpha_viewer",
    version="0.0.1",  
    description="alphafold viewer",
    long_description=long_description,
    long_description_content_type="text/markdown",  
    url="https://github.com/Intron7/alpha_viewer",  
    author="Severin Dicks",  # Optional
    author_email="severin.dicks@gmail.com",  # Optional
    package_dir={"":"alpha_viewer"},  # Optional
    packages=find_packages(where="alpha_viewer"),  # Required
    python_requires=">=3.8, <4",
    install_requires=required,  # Optional
    project_urls={  # Optional
        "Bug Reports": "https://github.com/Intron7/alpha_viewer/issues",
        "Source": "https://github.com/Intron7/alpha_viewer",
    },
)
