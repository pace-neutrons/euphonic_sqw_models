import versioneer
try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup

def get_euphonic_version():
    # gets the required euphonic version from `min_requirements.txt` file
    import os
    curdir = os.path.dirname(os.path.abspath(__file__))
    with open(os.path.join(curdir, 'min_requirements.txt')) as minreq:
        verstr = [req for req in minreq if 'euphonic' in req]
    return verstr[0].split('=')[1].strip()


setup(name='euphonic_horace',
      version=versioneer.get_version(),
      cmdclass=versioneer.get_cmdclass(),
      description='Phonon inelastic neutron spectra calculation using Horace and Euphonic',
      author='Rebecca Fair',
      author_email='rebecca.fair@stfc.ac.uk',
      url='https://github.com/pace-neutrons/euphonic_horace',
      packages=['euphonic_horace'],
      install_requires=['euphonic>=' + get_euphonic_version()],
     )
