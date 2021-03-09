import sys
import subprocess

pipcmd = [sys.executable, '-m', 'pip', 'install']

def clone_euphonic():
    subprocess.check_call(['git', 'clone', 'https://github.com/pace-neutrons/euphonic'])
    subprocess.check_call(pipcmd + ['.'], cwd='euphonic')


def install_requirements(requirements_file='min_requirements.txt'):
    with open(requirements_file, 'r') as req:
        packages = req.read().strip().split('\n')
    euphonic = [pp for pp in packages if 'euphonic' in pp]
    packages.pop(packages.index(euphonic[0]))
    try:
        subprocess.check_call(pipcmd + euphonic)
    except subprocess.CalledProcessError:
        clone_euphonic()
    subprocess.check_call(pipcmd + packages)


if __name__ == "__main__":
    install_requirements()
