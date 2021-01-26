import versioneer
try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup

euphonic_ver = '0.3.2'

setup(name='euphonic_horace',
      version=versioneer.get_version(),
      cmdclass=versioneer.get_cmdclass(),
      description='Phonon inelastic neutron spectra calculation using Horace and Euphonic',
      author='Rebecca Fair',
      author_email='rebecca.fair@stfc.ac.uk',
      url='https://github.com/pace-neutrons/euphonic_horace',
      packages=['euphonic_horace'],
      install_requires=['euphonic>=' + euphonic_ver],
     )
