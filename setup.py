"""Setup script for chimeric contig detector."""

from setuptools import setup, find_packages

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

with open("requirements.txt", "r", encoding="utf-8") as fh:
    requirements = [line.strip() for line in fh if line.strip() and not line.startswith("#")]

setup(
    name="chimeric-contig-detector",
    version="1.0.0",
    author="Megan Johnson",
    author_email="noreply@example.com",
    description="Detect chimeric contigs in viral metagenomic co-assemblies using read-pair analysis",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/megjohnson1999/chimera-break",
    packages=find_packages(),
    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
    ],
    python_requires=">=3.8",
    install_requires=requirements,
    entry_points={
        "console_scripts": [
            "detect-chimeras=detect_chimeras:main",
        ],
    },
    include_package_data=True,
    package_data={
        "chimeric_detector": ["example_config.yaml"],
    },
)