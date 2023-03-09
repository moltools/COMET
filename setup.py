import setuptools


setuptools.setup(
    name="COMET",
    version="0.0.1",
    author="",
    author_email="",
    description="",
    url="",
    install_requires=["rdkit>=2022.9.2", "scipy>=1.9.3"],
    packages=setuptools.find_packages(),
    python_requires=">=3.9",
    entry_points={"console_scripts": ["COMET = COMET.main:main"]}
)