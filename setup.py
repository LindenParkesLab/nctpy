import setuptools

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setuptools.setup(
    name="network-control",
    version="0.0.1",
    author="Jennifer Stiso",
    author_email="danibassettlab@gmail.com",
    description="Python implementation of concepts from network control theory",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/BassettLab/control_package",
    project_urls={
        "Bug Tracker": "hhttps://github.com/BassettLab/control_package/issues",
    },
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    package_dir={"": "src"},
    packages=setuptools.find_packages(where="src"),
    python_requires=">=3.6",
)