from setuptools import setup

config = {
    'include_package_data': True,
    'description': 'GenomeDISCO',
    'download_url': 'https://github.com/kundajelab/genomedisco',
    'version': '1.0.0',
    'packages': ['genomedisco','examples'],
    'setup_requires': [],
    'install_requires': ['numpy>=1.9', 'matplotlib>=1.5.0'],
    'scripts': [],
    'entry_points': {'console_scripts': ['genomedisco = genomedisco.__main__:main']},
    'name': 'genomedisco',
}

if __name__== '__main__':
    setup(**config)
