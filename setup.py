import setuptools


setuptools.setup(
    name="comet",
    version="0.0.1",
    author="",
    author_email="",
    description="",
    url="",
    install_requires=["rdkit>=2022.9.2", "scipy>=1.9.3", "tqdm"],
    package_dir={"": "src"},
    packages=["comet"],
    python_requires=">=3.9",
    entry_points={"console_scripts": ["comet = main:main"]}
)