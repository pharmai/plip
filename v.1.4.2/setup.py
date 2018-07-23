from setuptools import setup, find_packages, Extension

setup(
    name='plip',
    packages=['plip.module','plip'],
    version='1.4.2',
    description='PLIP - Fully automated protein-ligand interaction profiler',
    author='Sebastian Salentin',
    author_email='sebastian.salentin@tu-dresden.de',
    url='https://github.com/ssalentin/plip',
    keywords=['plip'],
    include_package_data=True,
    entry_points={'console_scripts': ['plip = plip.plipcmd:main']},
    zip_safe=False,
    classifiers=[
        '',
    ]
)
