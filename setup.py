from setuptools import setup, find_packages

version = {}
with open("src/_version.py") as version_file:
    exec(version_file.read(), version)

setup(
    name="buscoprot",
    version=version["__version__"],
    author="neng huang",
    license="Licensed under the Apache License. See LICENSE file.",
    author_email="neng@ds.dfci.harvard.edu",
    long_description="A BUSCO replacement based on miniprot to assessing genome assembly and annotation completeness "
    "with Benchmarking Universal Single-Copy Orthologs ",
    url="",
    download_url="",
    entry_points={
            'console_scripts': [
                'buscoprot=src.BuscoprotRunner:main',
                ],
            },
    platforms="Unix like",
    zip_safe=False,
    packages=["src"],
)
