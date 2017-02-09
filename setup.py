from setuptools import setup

config = {
    'include_package_data': True,
    'description': 'Genome Disco: DIfferences in Smoothed COntact maps',
    'download_url': 'https://github.com/kundajelab/genomedisco',
    'version': '0.1',
    'packages': ['genomedisco', 'genomedisco/comparison_types'],
    'setup_requires': [],
    'install_requires': ['numpy>=1.9', 'matplotlib<=1.5.3', 'scikit-learn'],
    'scripts': ['scripts/disco.py', 'scripts/annotate_baits.py', 'scripts/genomedisco_GenomewideIntraChromosomal.sh','scripts/genomedisco_multiChromosome.sh'],
    'name': 'genomedisco'
}

if __name__== '__main__':
    setup(**config)
