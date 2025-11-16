from setuptools import setup, find_packages
from pybind11.setup_helpers import Pybind11Extension, build_ext

# 精度優先で速度を狙うコンパイル／リンクオプション（fast-math は無効）
compile_flags = [
    "-O3",
    "-march=native",
    "-mtune=native",
    "-funroll-loops",
    "-fomit-frame-pointer",
    "-fno-fast-math",
    "-fvisibility=hidden",           # 不要シンボルを隠して最適化しやすくする
    "-fvisibility-inlines-hidden",
    "-fstrict-aliasing",             # ポインタ別名最適化を許可
]
link_flags = [
    "-flto",                         # LTO で関数間最適化
    "-Wl,-dead_strip",               # 未使用シンボルの除去（macOS系）
    "-Wl,-dead_strip_dylibs",
]
macros = [("NDEBUG", None)]  # アサート無効化で実行時オーバーヘッド低減

ext_modules = [
    Pybind11Extension(
        "triaccel._core",                   # ← サブモジュールとしてビルド
        ["triaccel/_core.cpp"],
        cxx_std=17,
        extra_compile_args=compile_flags,
        extra_link_args=link_flags,
        define_macros=macros,
    ),
]

setup(
    name="triaccel",
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
