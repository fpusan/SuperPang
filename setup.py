from setuptools import setup, Extension
import numpy
import setuptools

extensions = [Extension("superpang.lib.Compressor", sources=["src/superpang/lib/Compressor.pyx"]),
              Extension("superpang.lib.cutils", sources=["src/superpang/lib/cutils.pyx"]),
              Extension("superpang.lib.vtools", sources=["src/superpang/lib/vtools.pyx"])]

setup(
    setup_requires = ['cython', 'numpy', 'setuptools>=61.0'],
    ext_modules = extensions,
    version = open('src/superpang/VERSION').read().strip(),
    include_dirs = [numpy.get_include()],
    )


# pushing stuff to pip :
# rm -r src/build dist src/SuperPang.egg-info src/superpang/lib/*c
# python3 -m build
# python3 -m twine upload --repository pypi dist/*tar.gz
# rm -r src/build dist src/SuperPang.egg-info src/superpang/lib/*c
