import ez_setup
ez_setup.use_setuptools()

# Parse version string
import os
this_directory = os.path.dirname(os.path.abspath(__file__))
version_file = os.path.join(this_directory, "FastProject", "_version.py")
exec(open(version_file).read())

from setuptools import setup, find_packages

setup(
    name="FastProject",
    version=__version__,
    packages=find_packages(),

    entry_points={'console_scripts': ['fastproject = FastProject.CLI:entry']},

    install_requires=['numpy>=1.11', 'scipy>=0.18.0',
                      'scikit-learn>=0.18.1', 'pandas>=0.19.0'],

    include_package_data=True,

    author="David DeTomaso",
    author_email="David.DeTomaso@berkeley.edu",
    description="A python module for evaulating dimensionality reduction techniques on gene expression data.",
    keywords="Bioinformatics dimensionality reduction projection rna-seq",
    url="https://github.com/YosefLab/FastProject",
    license="Free for Academic and Not-For-Profit Use, See Repository for Details"
)
