from setuptools import setup, find_packages

version = {}
with open("_version.py") as version_file:
    exec(version_file.read(), version)

setup(
    name="compleasm",
    version=version["__version__"],
    author="Neng Huang",
    license="Licensed under the Apache License. See LICENSE file.",
    author_email="neng@ds.dfci.harvard.edu",
    long_description="Compleasm: a faster and more accurate reimplementation of BUSCO",
    url="https://github.com/huangnengCSU/compleasm",
    download_url='https://github.com/huangnengCSU/compleasm/archive/{}.tar.gz'.format(version),
    entry_points={
            'console_scripts': [
                'compleasm=compleasm:main',
                ],
            },
    platforms="Unix like",
    zip_safe=False,
    packages=["."],
)