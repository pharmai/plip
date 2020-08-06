from setuptools import setup, find_packages

from plip.basic import config

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
      install_requires=[
          'openbabel',
          'numpy',
          'lxml'
      ],
      entry_points={
          "console_scripts": [
              "plip = plip.plipcmd:main"
          ]
      },
      zip_safe=False)
