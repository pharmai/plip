from setuptools import setup, find_packages
from setuptools.command.build_ext import build_ext
from setuptools.command.install import install
from distutils.command.build import build

from plip.basic import config

class CustomBuild(build):
    """Ensure build_ext runs first in build command."""
    def run(self):
        self.run_command('build_ext')
        build.run(self)

class CustomInstall(install):
    """Ensure build_ext runs first in install command."""
    def run(self):
        self.run_command('build_ext')
        install.run(self)

class CustomBuildExt(build_ext):
    """ A workaround to build openbabel python3 bindings, for details see:
    https://github.com/openbabel/openbabel/issues/2408
    """

    def run(self):
        import requests
        import tarfile
        import tempfile
        import shutil
        import fileinput
        import subprocess
        import sys

        openbabel_pypi_url='https://files.pythonhosted.org/packages/9d/3f/f08f5d1422d74ed0e1e612870b343bfcc26313bdf9efec9165c3ea4b3ae2/openbabel-3.1.1.1.tar.gz'

        def install(package):
            subprocess.check_call([sys.executable, '-m', 'pip', 'install', '--user', package])

        print (f"Downloading openbabel package from : {openbabel_pypi_url}")
        obtar=requests.get(openbabel_pypi_url)
        obtmpdir = tempfile.mkdtemp()
        obtmp = obtmpdir+'/openbabel-3.1.1.1.tar'
        open(obtmp,'wb').write(obtar.content)
        print(f"Saving openbable tar.gz to {obtmpdir}")
        versfile = obtmpdir+'/openbabel-3.1.1.1/openbabel/__init__.py'

        with tarfile.open(obtmp,mode='r') as tf:
            tf.extractall(obtmpdir)

        print ('Fix versions: s/3.1.1.1/3.1.1/ to make StrictVersion() happy')
        print ('See https://github.com/openbabel/openbabel/issues/2408 for more details')
        with fileinput.input(files=versfile,inplace=True) as f:
            for line in f:
                op = line.replace('__version__ = "3.1.1.1"', '__version__ = "3.1.1"')
                print(op, end='')

        install(obtmpdir+'/openbabel-3.1.1.1/')
        print (f"Cleanup tmpdir: {obtmpdir}")
        shutil.rmtree(obtmpdir)

setup(name='plip',
        version=config.__version__,
        description='PLIP - Fully automated protein-ligand interaction profiler',
        classifiers=[
            'Development Status :: 5 - Production/Stable',
            'Intended Audience :: Science/Research',
            'Natural Language :: English',
            'License :: OSI Approved :: GNU General Public License v2 (GPLv2)',
            'Programming Language :: Python :: 3.6',
            'Topic :: Scientific/Engineering :: Bio-Informatics'
            ],
        url='https://github.com/pharmai/plip',
        author='PharmAI GmbH',
        author_email='hello@pharm.ai',
        license='GPLv2',
        packages=find_packages(),
        scripts=['plip/plipcmd.py'],
        cmdclass={'build': CustomBuild, 'build_ext': CustomBuildExt, 'install': CustomInstall},
        install_requires=[
            'requests',
            'numpy',
            'lxml'
            ],
        entry_points={
            "console_scripts": [
                "plip = plip.plipcmd:main"
                ]
            },
        zip_safe=False)
