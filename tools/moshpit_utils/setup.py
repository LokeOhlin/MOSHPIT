import setuptools
import subprocess
import os
from setuptools import setup, find_packages

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setuptools.setup(
    name="moshpit_utils",
    version="0.0.1",
    author="Loke Ohlin",
    author_email="loke.lonnblad@gmail.com",
    description="A package containing functions used for reducing the outputs used by moshpid",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://bitbucket.org/LokeLOhlin/moshpit",
    install_requires=["h5py", 
                      "numpy", 
                      "matplotlib"],
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    package_dir={"": "./"},
    packages=setuptools.find_packages(where="./"),
    python_requires=">=3.6",
)

# compile code bits
cwd = os.getcwd()
os.chdir("moshpit_utils/test_utils/DustyWave/")
subprocess.run("make")
os.chdir(cwd)


