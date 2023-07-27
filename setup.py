from Cython.Build import cythonize
import setuptools

extensions = [setuptools.Extension("superpang.lib.Compressor", sources=["src/superpang/lib/Compressor.pyx"]),
              setuptools.Extension("superpang.lib.cutils", sources=["src/superpang/lib/cutils.pyx"]),
              setuptools.Extension("superpang.lib.vtools", sources=["src/superpang/lib/vtools.pyx"])]

setuptools.setup(
    ext_modules = cythonize(extensions),
    setup_requires = ['cython'],
    version = open('src/superpang/VERSION').read().strip()
    )


# pushing stuff to pip :
# rm -r src/build dist src/SuperPang.egg-info src/superpang/lib/*c
# python3 -m build
# python3 -m twine upload --repository pypi dist/*tar.gz
# rm -r src/build dist src/SuperPang.egg-info src/superpang/lib/*c
