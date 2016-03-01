import ez_setup
ez_setup.use_setuptools()


from setuptools import setup, find_packages

setup(
	name = "FastProject",
	version = "1.0.2",
	packages = find_packages(),

	entry_points = { 'console_scripts': ['fastproject = FastProject.__main__:entry']},
	
	install_requires = ['numpy>=1.9','scipy>=0.14.0','scikit-learn>=0.15.2','pandas>=0.16.2'],
	
	include_package_data = True,
	
	author = "David DeTomaso",
	author_email = "David.DeTomaso@berkeley.edu",
	description = "A python module for evaulating dimensionality reduction techniques on gene expression data.",
	keywords = "Bioinformatics dimensionality reduction projection rna-seq",
	url = "https://github.com/YosefLab/FastProject",
	)
	
