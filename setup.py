"""
Protein-Ligand Interaction Profiler - Analyze and visualize protein-ligand interactions in PDB files.
setup.py - Setup configuration file for pip, etc.
Copyright 2014 Sebastian Salentin

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
"""

from setuptools import setup

setup(name='plip',
      version='1.2.2',
      description='PLIP - Fully automated protein-ligand interaction profiler',
      classifiers=[
          'Development Status :: 5 - Production/Stable',
          'Intended Audience :: Science/Research',
          'Natural Language :: English',
          'License :: OSI Approved :: Apache Software License',
          'Programming Language :: Python :: 2.7',
          'Topic :: Scientific/Engineering :: Bio-Informatics'
      ],
      url='https://github.com/ssalentin/plip',
      author='Sebastian Salentin',
      author_email='sebastian.salentin@biotec.tu-dresden.de',
      license='Apache',
      packages=['plip', 'plip/modules'],
      scripts=['plip/plipcmd'],
      install_requires=[
          'openbabel',
          'numpy',
          'lxml',
          'pymol',
      ],
      zip_safe=False)
