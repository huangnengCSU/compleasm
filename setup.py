from setuptools import setup, find_packages

version = {}
with open("_version.py") as version_file:
    exec(version_file.read(), version)

setup(
    name="minibusco",
    version=version["__version__"],
    author="Neng Huang",
    license="Licensed under the Apache License. See LICENSE file.",
    author_email="neng@ds.dfci.harvard.edu",
    long_description="miniBUSCO: a faster and more accurate reimplementation of BUSCO",
    url="https://github.com/huangnengCSU/minibusco",
    download_url='https://github.com/huangnengCSU/minibusco/archive/v{}.tar.gz'.format(version),
    entry_points={
            'console_scripts': [
                'minibusco=minibusco:main',
                ],
            },
    platforms="Unix like",
    zip_safe=False,
    packages=["."],
)