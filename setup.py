import setuptools

with open("README.rst", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setuptools.setup(
    name="nctpy",
    version="1.0.0.post0",

    author="Linden Parkes, Jason Kim, and Jennifer Stiso",
    author_email="lindenparkes@gmail.com, jinsu1@seas.upenn.edu, jeni.stiso@gmail.com",

    description="Python implementation of concepts from network control theory",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/BassettLab/nctpy",
    project_urls={
        "Bug Tracker": "https://github.com/BassettLab/nctpy/issues",
    },
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    package_dir={"": "src"},
    packages=setuptools.find_packages(where="src"),
    python_requires=">=3.9",
    install_requires=['numpy', 'scipy', 'tqdm', 'statsmodels', 'seaborn', 'nibabel', 'nilearn']
)
