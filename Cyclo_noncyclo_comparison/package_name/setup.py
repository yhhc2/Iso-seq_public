from setuptools import setup, find_packages

setup(
    name="package_name",
    version="0.1",
    packages=find_packages(),
    install_requires=[
        "pandas",
        "numpy",
        "matplotlib",
        "statsmodels",
        "seaborn",
        "scipy",
    ],
    author="Hank Cheng",
    description="Package for calculating test statistics and z-scores for isoforms.",
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
    ],
)
