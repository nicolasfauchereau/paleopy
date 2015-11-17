#!/usr/bin/env python

from setuptools import setup, Command

from pip.req import parse_requirements

# parse_requirements() returns generator of pip.req.InstallRequirement objects
install_reqs = parse_requirements('requirements.txt')

# reqs is a list of requirement
# e.g. ['django==1.5.1', 'mezzanine==1.4.6']
reqs = [str(ir.req) for ir in install_reqs]

# TBD fix dependencies
install_requires = reqs
extras_requires = {
	'tests': ['pytest']
}
# TBD get version from module __init__
version = '0.1'
packages = ['paleopy']

class PyTest(Command):
    user_options = []
    def initialize_options(self):
        pass

    def finalize_options(self):
        pass

    def run(self):
        import subprocess
        import sys
        errno = subprocess.call([sys.executable, 'runtests.py'])
        raise SystemExit(errno)

setup(
    name = "paleopy",
    description = 'paleopy',
    long_description = 'paleopy',
    version = version,
    url = 'https://github.com/nicolasfauchereau/paleopy',
    license = 'All rights reserved',
    platforms = ['unit', 'linux', 'osx', 'cygwin', 'win32'],
    author = 'Nicolas Fauchereau',
    author_email = 'nicolas dot fauchereau at niwa co nz',
    entry_points = {},
    install_requires = install_requires,
    extras_require = extras_requires,
    packages = packages,
    zip_safe = False,
    cmdclass = {'test': PyTest}
)