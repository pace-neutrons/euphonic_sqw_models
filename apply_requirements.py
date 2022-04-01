"""
This script installs test dependencies, and the minimum required
version of Euphonic via pip for use in testing. If a breaking change
has been applied to Euphonic, the min requirements for Euphonic should
be something like euphonic>0.6.0 (note: no '='), in which case the
latest specified branch of Euphonic will be cloned, installed and
tested against
"""

import sys
import os
import subprocess
from argparse import ArgumentParser

pipcmd = [sys.executable, '-m', 'pip', 'install']

def clone_euphonic(branch='master'):
    subprocess.check_call(
        ['git', 'clone', 'https://github.com/pace-neutrons/Euphonic'])
    subprocess.check_call(['git', 'checkout', branch], cwd='Euphonic')


def install_requirements(requirements_file=None,
                         extras=[], branch='master'):
    if requirements_file is None:
        requirements_file = os.path.join(os.path.dirname(__file__),
                                         'min_requirements.txt')
    with open(requirements_file, 'r') as req:
        packages = req.read().strip().split('\n')
    euphonic = [pp for pp in packages if 'euphonic' in pp][0]
    packages.pop(packages.index(euphonic))
    euphonic_ver = f'>{euphonic.split(">")[-1]}'
    euphonic_extras = f'[{",".join(extras)}]'
    try:
        subprocess.check_call(
            pipcmd + ['euphonic' + euphonic_extras + euphonic_ver])
    except subprocess.CalledProcessError:
        clone_euphonic(branch=branch)
        subprocess.check_call(pipcmd + ['.' + euphonic_extras],
                              cwd='Euphonic')
    subprocess.check_call(pipcmd + packages)


if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument(
        '--branch', default='master',
        help='The Euphonic branch to clone')
    parser.add_argument(
        '--extras', nargs='*', default=[],
        help='The Euphonic "extras" to install e.g. phonopy_reader')
    args = parser.parse_args()

    install_requirements(branch=args.branch, extras=args.extras)
