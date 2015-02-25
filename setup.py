import ez_setup
ez_setup.use_setuptools()


from setuptools import setup, find_packages

setup(
	name = "FastProject",
	version = "0.7.1",
	packages = find_packages(),
	
	install_requires = ['numpy>=1.9','matplotlib>=1.4.0','scikit-learn>=0.15.2'],
	
	include_package_data = True,
	
	author = "David DeTomaso",
	author_email = "David.DeTomaso@berkeley.edu",
	description = "A python module for evaulating dimensionality reduction techniques on gene expression data.",
	keywords = "Bioinformatics dimensionality reduction projection",
	url = "https://github.com/deto/FastProject",
	)
	
