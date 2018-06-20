from setuptools import find_packages, setup


setup(
    name="pylprotpredictor",
    version="0.1.0",
    author="Berenice Batut",
    author_email="berenice.batut@gmail.com",
    description=("A tool to predict CDS including Pyrrolise-using CDS"),
    license="Apache-2.0",
    keywords="",
    url="https://github.com/bebatut/PylProtPredictor ",
    packages=find_packages(),
    scripts=['bin/pylprotpredictor'],
    classifiers=[
        "Development Status :: 2 - Pre-Alpha",
        "License :: OSI Approved :: Apache-2.0",
        "Programming Language :: Python :: 3",
        "Intended Audience :: Science/Research",
        "Operating System :: OS Independent",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.4",
        "Programming Language :: Python :: 3.5",
        "Programming Language :: Python :: 3.6"
    ],
    extras_require={
        'testing': ["pytest"],
    },
    setup_requires=['pytest-runner'],
    tests_require=['pytest'],
    install_requires=[
        'biopython',
        'flake8',
        'codecov',
        'pandas',
        'pytest-cov',
        'requests',
        'Sphinx',
        'sphinx_rtd_theme']
)