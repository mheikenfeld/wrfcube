from setuptools import setup

setup(name='wrfcube',
      use_scm_version=True,
      setup_requires=['setuptools_scm'],
      description='Load WRF output into iris cubes',
      url='http://github.com/mheikenfeld/wrfcube',
      author='Max Heikenfeld',
      author_email='maxheikenfeld@web.de',
      license='BSD-3-Clause',
      packages=['wrfcube'],
      install_requires=[],
      zip_safe=False)
