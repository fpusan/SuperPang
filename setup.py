import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="SuperPang",
    version=open('VERSION').read().strip(),
    author="Fernando Puente-Sánchez",
    author_email="fernando.puente.sanchez@slu.se",
    description="Non-redundant pangenome assemblies from multiple genomes or bins",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/fpusan/SuperPang",
    # BSD 3-Clause License:
    # - http://choosealicense.com/licenses/bsd-3-clause
    # - http://opensource.org/licenses/BSD-3-Clause
    license='BSD',
    packages= setuptools.find_packages(),
    install_requires = [
    "numpy", "mOTUlizer==0.2.4", "mappy", # and graph-tool, but we need conda for that
    ],
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: BSD License",
        "Operating System :: OS Independent",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics"
    ],
    keywords='bioinformatics assembly metagenomics microbial-genomics genomics',
    python_requires='>=3.8',
    scripts = ['SuperPang/scripts/SuperPang.py','SuperPang/scripts/homogenize.py', 'SuperPang/scripts/condense-edges.py', 'SuperPang/scripts/test-SuperPang.py'],
    #package_data = {'SuperPang': ['SuperPang/test_data/*']},
    include_package_data = True,
)

# pushing stuff to pip : python3 setup.py sdist bdist_wheel
# python3 -m twine upload --repository pypi dist/*
