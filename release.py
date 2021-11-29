import argparse
import json
import os
import re
import requests
import warnings

from packaging import version

from euphonic_sqw_models import __version__


def main():
    parser = get_parser()
    args = parser.parse_args()

    print(args)
    test = not args.notest
    if args.github:
        release_github(test)

def check_euphonic_version():
    """
    Check euphonic_sqw_models has a release version of Euphonic as a
    requirement, and raise an Exception if not. If so, check that this
    version is Euphonic's latest release version and raise a warning if
    not
    """
    req_file = 'min_requirements.txt'
    with open(req_file) as f:
        lines = f.readlines()
    found_euphonic_requirement = False
    for line in lines:
        if 'euphonic' in line:
            found_euphonic_requirement = True
            if not '=' in line:
                raise Exception(f'Release versions of euphonic_sqw_models '
                                f'should be compatible with a Euphonic '
                                f'release version! Expected {line.strip()} to '
                                f'contain a "="')
            euphonic_ver = line.strip()[-5:]
            response = requests.get(
                    'https://pypi.org/pypi/euphonic/json')
            euphonic_rel_vers = [version.parse(x)
                                 for x in list(response.json()['releases'].keys())]
            latest_euphonic_release = max(euphonic_rel_vers)
            if str(latest_euphonic_release) != euphonic_ver:
                warnings.warn(f'The Euphonic version specified in {req_file} '
                              f'({euphonic_ver}) is not the latest release '
                              f'version ({latest_euphonic_release}). Are you '
                              f'sure this is correct?')
    if not found_euphonic_requirement:
        raise Exception(f'Euphonic requirement not found in {req_file}')


def release_github(test=True):
    check_euphonic_version()
    with open('CHANGELOG.rst') as f:
        changelog = f.read()
    eu_sqw_models_ver = 'v' + __version__
    changelog_ver = re.findall('\n`(v\d+\.\d+\.\S+)\s', changelog)[0]
    if eu_sqw_models_ver != changelog_ver:
        raise Exception((
            f'VERSION and CHANGELOG.rst version mismatch!\n'
            f'VERSION: {eu_sqw_models_ver}\nCHANGELOG.rst: '
            f'{changelog_ver}'))
    desc = re.search('`v\d+\.\d+\.\S+.*?^-+\n(.*?)^`v', changelog,
                     re.DOTALL | re.MULTILINE).groups()[0].strip()

    payload = {
        "tag_name": changelog_ver,
        "target_commitish": "main",
        "name": changelog_ver,
        "body": desc,
        "draft": False,
        "prerelease": False
    }
    if test:
        print(payload)
    else:
        response = requests.post(
            'https://api.github.com/repos/pace-neutrons/euphonic_sqw_models/releases',
            data=json.dumps(payload),
            headers={"Authorization": "token " + os.environ["GITHUB_TOKEN"]})
        print(response.text)

def get_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '--github',
        action='store_true',
        help='Release on Github')
    parser.add_argument(
        '--notest',
        action='store_true',
        help='Actually send/upload')
    return parser


if __name__ == '__main__':
    main()

