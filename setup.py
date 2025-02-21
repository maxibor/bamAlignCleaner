from setuptools import find_packages, setup

from bamaligncleaner import __version__

setup(
    name="bamAlignCleaner",
    version=__version__,
    description="Removes unaligned references in BAM alignment file",
    long_description=open("README.md").read(),
    url="https://github.com/maxibor/bamAlignCleaner",
    long_description_content_type="text/markdown",
    license="MIT",
    python_requires=">=3.6",
    install_requires=["pysam", "tqdm", "click", "numpy"],
    packages=find_packages(include=["bamaligncleaner"]),
    entry_points={"console_scripts": ["bamAlignCleaner= bamaligncleaner.cli:cli"]},
)
