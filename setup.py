from setuptools import setup, find_packages

setup(
    name='TargetSiteDuplicationFinder',
    version='1.0.0',
    description='Script for analysing alignment files generated from BWA-MEM2 for transposable element insertions',
    author='Andreas Bachler',
    author_email='Andy.Bachler@gmail.com',
    url='https://github.com/Andy-B-123/TargetSiteDuplicationFinder',
    packages=find_packages(),
    install_requires=[
        'tqdm',
        'pysam',
        'numpy',
        'icecream',
        'pandas'
        
    ],
    entry_points={
        'console_scripts': [
            'your_script_name = TSD_parser.cluster_wrapped:main',
        ],
    },
)
