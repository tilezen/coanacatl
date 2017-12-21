from setuptools import setup
from distutils.extension import Extension
import os.path

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
    ext_modules=[extension],
    test_suite='test',
)
