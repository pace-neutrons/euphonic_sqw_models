import sys
import subprocess

pipcmd = [sys.executable, '-m', 'pip', 'install']

def clone_euphonic(suffix=''):
    subprocess.check_call(['git', 'clone', 'https://github.com/pace-neutrons/euphonic'])
    subprocess.check_call(pipcmd + ['.' + suffix], cwd='euphonic')


def install_requirements(requirements_file='min_requirements.txt'):
    with open(requirements_file, 'r') as req:
        packages = req.read().strip().split('\n')
    euphonic = [pp for pp in packages if 'euphonic' in pp][0]
    packages.pop(packages.index(euphonic))
    if '[' in euphonic and ']' in euphonic:
        suffix = f"[{euphonic.split('[')[1].split(']')[0]}]"
    else:
        suffix = ''
    try:
        subprocess.check_call(pipcmd + [euphonic])
    except subprocess.CalledProcessError:
        clone_euphonic(suffix)
    subprocess.check_call(pipcmd + packages)


if __name__ == "__main__":
    install_requirements()
