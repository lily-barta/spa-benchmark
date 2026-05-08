from setuptools import setup, find_packages

setup(
    name="spa-benchmark",
    version="0.1",
    packages=find_packages(where="src"),
    package_dir={"": "src"},
    install_requires=[
        "pandas",
        "matplotlib", 
        "qulacs",
        "pyscf",
        "openfermion",
        "project-sunrise==0.2"
    ],
    url="https://github.com/lily-barta/spa-benchmark",
    author="Lily Barta",
)
