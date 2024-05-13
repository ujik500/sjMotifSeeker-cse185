from setuptools import *

setup(
    name = "sjMotifSeeker",
    version = "1.0",
    description = "CSE 185 Project: Motif Finding",
    author = "Jack Kissinger & Sujana Sreenivasan",
    author_email = "jkissinger@ucsd.edu",
    packages = find_packages(),
    entry_points = {
        "console_scripts": ["sjMotifSeeker=source.sjMotifSeeker:main"],
    },
)