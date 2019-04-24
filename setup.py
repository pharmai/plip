"""
Protein-Ligand Interaction Profiler - Analyze and visualize protein-ligand interactions in PDB files.
setup.py - Setup configuration file for pip, etc.
"""

from setuptools import setup

setup(name='plip',
      version='1.4.5',
      description='PLIP - Fully automated protein-ligand interaction profiler',
      classifiers=[
          'Development Status :: 5 - Production/Stable',
          'Intended Audience :: Science/Research',
          'Natural Language :: English',
          'License :: OSI Approved :: GNU General Public License v2 (GPLv2)',
          'Programming Language :: Python :: 2.7',
          'Programming Language :: Python :: 3.5',
          'Topic :: Scientific/Engineering :: Bio-Informatics'
      ],
      url='https://github.com/ssalentin/plip',
      author='Sebastian Salentin, Joachim Haupt',
      author_email='joachim.haupt@tu-dresden.de',
      license='GPLv2',
      packages=['plip', 'plip/modules'],
      scripts=['plip/plipcmd.py'],
      install_requires=[
          'openbabel',
          'numpy',
          'lxml',
          'future'
      ],
      zip_safe=False)
