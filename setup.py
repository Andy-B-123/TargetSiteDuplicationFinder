from setuptools import setup, find_packages

setup(
    name='your_script_name',
    version='1.0.0',
    description='Description of your script',
    author='Your Name',
    author_email='your@email.com',
    url='https://github.com/yourusername/yourrepository',
    packages=find_packages(),
    install_requires=[
        'tqdm',
        'pysam',
        'numpy'
    ],
    entry_points={
        'console_scripts': [
            'your_script_name = your_module_name:main',
        ],
    },
)
