from setuptools import setup, find_packages
from pybind11.setup_helpers import Pybind11Extension, build_ext

ext_modules = [
    Pybind11Extension(
        "triaccel._core",                   # ← サブモジュールとしてビルド
        ["triaccel/_core.cpp"],
        cxx_std=17,
        extra_compile_args=["-O3", "-ffast-math"],
    ),
]

setup(
    name="triaccel",
    version="0.4.0",
    packages=find_packages(),
    ext_modules=ext_modules,
    cmdclass={"build_ext": build_ext},
    include_package_data=True,
    package_data={
        "triaccel": [
            "py.typed",
            "*.pyi",
            "detail/*.hpp",
        ]
    },
    zip_safe=False,
)
