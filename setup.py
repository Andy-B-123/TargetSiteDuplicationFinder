from setuptools import setup, find_packages

setup(
    name='your_script_name',
    version='1.0.0',
    description='Description of your script',
    author='Your Name',
    author_email='your@email.com',
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
