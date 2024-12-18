# setup.py

import os
import platform
import subprocess
import sys
from distutils.version import LooseVersion

from setuptools import Extension, find_packages, setup
from setuptools.command.build_ext import build_ext


class CMakeExtension(Extension):
    """
    A setuptools Extension for CMake projects.
    """

    def __init__(self, name, sourcedir=''):
        super().__init__(name, sources=[])
        self.sourcedir = os.path.abspath(sourcedir)


class CMakeBuild(build_ext):
    """
    A custom build extension for building CMake projects within setuptools.
    """

    def run(self):
        # Ensure CMake is installed
        try:
            subprocess.check_output(['cmake', '--version'])
        except OSError:
            raise RuntimeError("CMake must be installed to build the extensions.")

        # Build each CMakeExtension
        for ext in self.extensions:
            self.build_extension(ext)

    def build_extension(self, ext):
        # Absolute path to the directory where the extension will be placed
        # extdir = os.path.abspath(os.path.dirname(self.get_ext_fullpath(ext.name)))

        # Create the build temporary directory if it doesn't exist
        if not os.path.exists(self.build_temp):
            os.makedirs(self.build_temp)

        # CMake configuration arguments
        cmake_args = [f"-DCMAKE_BUILD_TYPE={'Debug' if self.debug else 'Release'}", f"-DPYTHON_EXECUTABLE={sys.executable}"]

        # CMake build arguments
        build_args = []
        cfg = 'Debug' if self.debug else 'Release'

        if platform.system() == "Windows":
            cmake_args += [f"-DCMAKE_LIBRARY_OUTPUT_DIRECTORY_{cfg.upper()}={self.build_temp}"]
            build_args += ['--', '/m']
        else:
            build_args += ['--', '-j2']  # Adjust '-j2' based on your CPU cores

        # Environment variables for the build
        env = os.environ.copy()
        env['CXXFLAGS'] = f'{env.get("CXXFLAGS", "")} -DVERSION_INFO="{self.distribution.get_version()}"'

        # Configure the project with CMake
        subprocess.check_call(['cmake', ext.sourcedir] + cmake_args, cwd=self.build_temp, env=env)

        # Build the project
        subprocess.check_call(['cmake', '--build', '.'] + build_args, cwd=self.build_temp)

        # Install the library (places EOBLib.so into pySEOBNREPHM/lib)
        subprocess.check_call(['cmake', '--install', '.'], cwd=self.build_temp)


def read_requirements():
    """
    Reads the requirements from requirements.txt.

    Returns:
        list: A list of dependencies.
    """
    here = os.path.abspath(os.path.dirname(__file__))
    with open(os.path.join(here, 'requirements.txt'), 'r') as f:
        return [line.strip() for line in f if line.strip() and not line.startswith('#')]


def read_requirements_dev():
    """
    Reads the development requirements from requirements-dev.txt.

    Returns:
        list: A list of development dependencies.
    """
    here = os.path.abspath(os.path.dirname(__file__))
    with open(os.path.join(here, 'requirements-dev.txt'), 'r') as f:
        return [line.strip() for line in f if line.strip() and not line.startswith('#')]


setup(
    name='pySEOBNREPHM',
    version='1.0.0',
    author='Your Name',
    author_email='your.email@example.com',
    description='pySEOBNREPHM Python Package with C Backend',
    url='https://github.com/yourusername/pySEOBNREPHM',
    packages=find_packages(),
    ext_modules=[CMakeExtension('pySEOBNREPHM')],
    cmdclass={
        'build_ext': CMakeBuild,
    },
    zip_safe=False,
    include_package_data=True,  # Include files specified in MANIFEST.in
    package_data={
        'pySEOBNREPHM': ['lib/libEOB.so'],
    },
    classifiers=[
        'Programming Language :: Python :: 3',
        'License :: OSI Approved :: MIT License',
        'Operating System :: OS Independent',
    ],
    python_requires='>=3.6',
    install_requires=read_requirements(),
    extras_require={
        'dev': read_requirements_dev(),
    },
)
