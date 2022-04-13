import setuptools
from setuptools.command.install import install
from setuptools.command.sdist import sdist
from sys import argv

class InstallCommand(install): # so that setup.py install will accept the "--cythonize" flag
    user_options = install.user_options + [
        ('cythonize', None, None),
    ]
    
    def initialize_options(self):
        install.initialize_options(self)
        self.cythonize = None

    def finalize_options(self):
        install.finalize_options(self)

    def run(self):
        install.run(self)

class SDistCommand(sdist):
    user_options = sdist.user_options + [
        ('cythonize', None, None),
    ]

    def initialize_options(self):
        sdist.initialize_options(self)
        self.cythonize = None

    def finalize_options(self):
        sdist.finalize_options(self)

    def run(self):
        sdist.run(self)


USE_CYTHON = "--cythonize" in argv # the --cythonize flag will make setup.py freak out unless we use it inside "sdist" or "install"

ext = '.pyx' if USE_CYTHON else '.c' # redistribute the c extension when uploading to pypi

extensions = [setuptools.Extension("SuperPang.lib.Compressor", sources=["SuperPang/lib/Compressor"+ext])]

if USE_CYTHON:
    from Cython.Build import cythonize
    extensions = cythonize(extensions)

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    cmdclass={'install': InstallCommand, 'sdist': SDistCommand},
    name="SuperPang",
    version=open('SuperPang/VERSION').read().strip(),
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
    ],                                    #  ... and minimus2 but that's not a python package
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
    include_package_data = True,
    ext_modules = extensions
)

# pushing stuff to pip :
# rm -r dist/* SuperPang.egg-info/ rm SuperPang/lib/*c SuperPang/lib/*so
# python3 setup.py sdist --cythonize
# python3 -m twine upload --repository pypi dist/*
# rm -r dist/* SuperPang.egg-info/ SuperPang/lib/*c
