from setuptools import setup

config = {
    'include_package_data': True,
    'description': 'GenomeDISCO',
    'download_url': 'https://github.com/kundajelab/genomedisco',
    'version': '1.0.0',
    'packages': ['genomedisco'],
    'setup_requires': [],
    'install_requires': ['numpy>=1.9', 'matplotlib>=1.5.0'],
    #'dependency_links': ["https://github.com/kundajelab/deeplift/tarball/v0.5.1-theano#egg=deeplift-0.5.1-theano",
    #                     "https://github.com/kundajelab/simdna/tarball/0.3#egg=simdna-0.3"],
    'scripts': [],
    'entry_points': {'console_scripts': ['genomedisco = genomedisco.__main__:main']},
    'name': 'genomedisco'
}

if __name__== '__main__':
    setup(**config)
