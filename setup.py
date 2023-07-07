from Cython.Build import cythonize
import setuptools

extensions = [setuptools.Extension("superpang.lib.Compressor", sources=["src/superpang/lib/Compressor.pyx"]),
              setuptools.Extension("superpang.lib.cutils", sources=["src/superpang/lib/cutils.pyx"]),
              setuptools.Extension("superpang.lib.vtools", sources=["src/superpang/lib/vtools.pyx"])]

setuptools.setup(
    ext_modules = cythonize(extensions),
    setup_requires = ['cython']
    )


# pushing stuff to pip :
# rm -r build/* dist/* src/SuperPang.egg-info/ src/superpang/lib/superPang/lib/*c
# python3 -m build
# python3 -m twine upload --repository pypi dist/*tar.gz
# rm -r build dist src/SuperPang.egg-info src/superPang/lib/*c
