from setuptools import setup

setup(name='wrfcube',
      version='0.1',
      description='Load WRF output into iris cubes',
      url='http://github.com/mheikenfeld/wrfcube',
      author='Max Heikenfeld',
      author_email='max.heikenfeld@physics.ox.ac.uk',
      license='GNU',
      packages=['wrfcube'],
      install_requires=['iris','numpy','netCDF4'],
      zip_safe=False)
