from setuptools import setup, find_packages

setup(
    name='infernommatic',  # This is the name your users will use for 'import' and 'pip install'
    version='0.0.1',
    packages=find_packages(), # Automatically finds your package directories
    install_requires=[
        'pandas>=2.2.2',
        'biopython>=1.84'
    ],
    author='Fernando Rossi',
    author_email='',
    description='This script was developed to address a problem at the SetuLab (University of Paulo), in which reads were still contaminated with adapters even after running multiple read trimming tools.',
    long_description=open('README.md').read(),
    long_description_content_type='text/markdown',
    url='https://github.com/FefoRossi/Infernommatic', # Your GitHub repo URL
    classifiers=[
        'Programming Language :: Python :: 3',
        'License :: OSI Approved :: MIT License',
        'Operating System :: OS Independent',
    ],
    python_requires='>=3.12.2', # Specify minimum Python version
    entry_points={'console_scripts': ['infernommatic = infernommatic.cli']}
)

