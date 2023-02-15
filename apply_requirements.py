"""
This script installs test dependencies, and the latest version of
Euphonic via pip for use in testing. If a breaking change
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
    # Only clone the first time
    if not os.path.exists('Euphonic'):
        subprocess.check_call(
            ['git', 'clone', 'https://github.com/pace-neutrons/Euphonic'])
    subprocess.check_call(['git', 'checkout', branch], cwd='Euphonic')


def install_requirements(requirements_file=None, extras=None,
                         branch=None, version=None):
    if branch is not None and version is not None:
        raise ValueError('Only one of branch or version should be provided')
    if requirements_file is None:
        requirements_file = os.path.join(os.path.dirname(__file__),
                                         'min_requirements.txt')
    with open(requirements_file, 'r') as req:
        packages = req.read().strip().split('\n')
    euphonic = [pp for pp in packages if 'euphonic' in pp][0]
    packages.pop(packages.index(euphonic))
    if '>=' not in euphonic and branch is None:
        # This means we're between versions, or a branch has been provided, clone Euphonic
        branch = 'master'
    if extras is None:
        ex_str = ''
    else:
        ex_str = f'[{",".join(extras)}]'
    if branch is not None:
        clone_euphonic(branch=branch)
        subprocess.check_call(pipcmd + [f'.{ex_str}'], cwd='Euphonic')
    else:
        if version is None:
            version_str = ''
        else:
            version_str = f'=={version}'
        subprocess.check_call(
            pipcmd + ['--upgrade'] + [f'euphonic{ex_str}{version_str}'])
    subprocess.check_call(pipcmd + packages)


if __name__ == "__main__":
    parser = ArgumentParser()
    euphonic_ver_group = parser.add_mutually_exclusive_group()
    euphonic_ver_group.add_argument(
        '--version', default=None,
        help='The Euphonic version to install, if not provided '
             'automatically uses most recent release')
    euphonic_ver_group.add_argument(
        '--branch', default=None,
        help='The Euphonic branch to clone, default master')
    parser.add_argument(
        '--extras', nargs='*', default=None,
        help='The Euphonic "extras" to install e.g. phonopy_reader')
    args = parser.parse_args()

    install_requirements(branch=args.branch, extras=args.extras,
                         version=args.version)
