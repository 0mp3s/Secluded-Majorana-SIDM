from setuptools import setup, find_packages

setup(
    name="secluded-majorana-sidm",
    version="1.0.0",
    description=(
        "Self-Interacting Dark Matter via a Secluded Scalar Mediator: "
        "Majorana Fermion with Velocity-Dependent Cross Sections"
    ),
    author="Omer P.",
    url="https://github.com/0mp3s/Secluded-Majorana-SIDM",
    license="MIT",
    packages=find_packages(),
    python_requires=">=3.10",
    install_requires=[
        "numpy",
        "scipy",
        "numba>=0.60",
        "matplotlib",
        "emcee>=3.1",
        "corner>=2.2",
        "pandas",
    ],
)
