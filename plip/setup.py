from setuptools import setup


setup(
    name='Protein-LigandInteractionProfiler',
    version='1.0.0',
    author='Sebastian Salentin',
    author_email='sebastian.salentin@biotec.tu-dresden.de',
    packages=['plip', 'plip.modules', 'docs'],
    scripts=['plip/plip-cmd.py'],
    url="",  # @todo Add url from github
    license='docs/LICENSE.txt',
    description='Characterize protein-ligand interactions from PDB files.',
    long_description=open('docs/README.rst').read(),
    install_requires=[
        "pymol",
        "openbabel",
        "numpy",
        "lxml"
    ],
)
