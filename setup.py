import versioneer
try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup

def get_required_versions():
    # gets required module versions from `min_requirements.txt` file
    import os
    curdir = os.path.dirname(os.path.abspath(__file__))
    with open(os.path.join(curdir, 'min_requirements.txt')) as minreq:
        reqs = minreq.read().splitlines()
    # Don't have to account for special case of Euphonic intermediate
    # versions here (e.g. euphonic>0.6.0) because that indicates a
    # dev/test version in which case apply_requirements.py should be
    # used. So just return reqs.
    return reqs


setup(name='euphonic_sqw_models',
      version=versioneer.get_version(),
      cmdclass=versioneer.get_cmdclass(),
      description='Phonon inelastic neutron spectra calculation using Horace and Euphonic',
      author='Rebecca Fair',
      author_email='rebecca.fair@stfc.ac.uk',
      url='https://github.com/pace-neutrons/euphonic_sqw_models',
      packages=['euphonic_sqw_models'],
      install_requires=get_required_versions(),
     )
