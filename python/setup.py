from setuptools import setup, find_packages

with open("README.md", "r") as fh:
    long_description = fh.read()

setup(
    name="filtering-and-smoothing",
    version="0.0.1",
    author="Adrien Corenflos",
    author_email="adrien.corenflos@aalto.fi",
    description="Companion Code for Bayesian Filtering and Smoothing book",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/EEA-sensors/Bayesian-Filtering-and-Smoothing",
    packages=find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    install_requires=["matplotlib>=3.2,<3.3",
                      "numpy>=1.19",
                      "scipy>=1.5"],
    python_requires='>=3.6',
)