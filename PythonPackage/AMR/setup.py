from setuptools import setup, find_packages

setup(
    name='AMR',
    version='2.1.1.9195',
    packages=find_packages(),
    install_requires=[
        'rpy2',
        'numpy',
        'pandas',
    ],
    author='Matthijs Berends',
    author_email='m.s.berends@umcg.nl',
    description='A Python wrapper for the AMR R package',
    long_description=open('README.md').read(),
    long_description_content_type='text/markdown',
    url='https://github.com/msberends/AMR',
    project_urls={
        'Bug Tracker': 'https://github.com/msberends/AMR/issues',
    },
    license='GPL 2',
    classifiers=[
        'Programming Language :: Python :: 3',
        'Operating System :: OS Independent',
    ],
    python_requires='>=3.6',
)
