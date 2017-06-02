from setuptools import setup

config = {
    'include_package_data': True,
    'description': 'Genome Disco: DIfferences in Smoothed COntact maps',
    'download_url': 'https://github.com/kundajelab/genomedisco',
    'version': '0.1.0',
    'packages': ['genomedisco', 'genomedisco/comparison_types'],
    'setup_requires': [],
    'entry_points': {'console_scripts': ['genomedisco = genomedisco.__main__:main']},
    'install_requires': ['numpy>=1.9', 'matplotlib<=1.5.3', 'scikit-learn'],
    'scripts': [],
    'name': 'genomedisco'
}

if __name__== '__main__':
    setup(**config)
