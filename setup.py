from setuptools import setup
from distutils.extension import Extension
import os.path

version_path = os.path.join(os.path.dirname(__file__), 'VERSION')
with open(version_path) as fh:
    version = fh.read().strip()

extension = Extension(
    "coanacatl",
    ["coanacatl/coanacatl.cpp"],
    include_dirs=[
        'vendor/mapbox/geometry.hpp/include',
        'vendor/mapbox/wagyu/include',
        'vendor/mapbox/protozero/include',
        'vendor/mapbox/vtzero/include',
    ],
    libraries=['boost_python', 'geos_c'],
    extra_compile_args=['-std=c++11'],
    language='c++11',
)

setup(
    name="coanacatl",
    version=version,
    ext_modules=[extension],
    test_suite='test',
)
