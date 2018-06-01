from setuptools import setup
setup(
	name='pyDNM',
	version='0.0.2.1',
	url='https://github.com/dantaki/pydnm',
	author='Danny Antaki, Aojie Lian',
	author_email='dantaki@ucsd.edu',
	packages=['pyDNM'],
	package_dir={'pyDNM': 'pyDNM/'},
	include_package_data=True,
	scripts= ['pyDNM/pydnm']
)