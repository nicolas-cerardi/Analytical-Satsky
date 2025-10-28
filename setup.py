from setuptools import setup, find_packages

setup(
   name='Analytical-Satsky',
   version='0.0.1',
   description='Analytical model of satellite exposure for telescopes',
   license="MIT",
   author='Nicolas Cerardi',
   author_email='nicolas.cerardi@epfl.ch',
   packages=find_packages(),  #same as name
   install_requires=['astropy'],
)