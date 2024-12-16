# setup.py

from setuptools import setup, find_packages
import os
import sys
from setuptools.command.install import install as _install

class InstallWithEOBLib(_install):
    """Custom install command to ensure EOBLib.so is included."""

    def run(self):
        # Run the standard install
        _install.run(self)
        # Ensure that EOBLib.so is copied to the package's lib directory
        # This is already handled via package_data, so additional steps may not be necessary
        # However, if you have post-install steps, add them here
        pass

setup(
    name='pySEOBNREPHM',
    version='1.0.0',
    author='Your Name',
    author_email='your.email@example.com',
    description='pySEOBNREPHM Python Package with C Backend',
    url='https://github.com/yourusername/pySEOBNREPHM',
    packages=find_packages(),
    include_package_data=True,  # Include files specified in MANIFEST.in
    package_data={
        'pySEOBNREPHM': ['lib/EOBLib.so'],
    },
    classifiers=[
        'Programming Language :: Python :: 3',
        'License :: OSI Approved :: MIT License',
        'Operating System :: OS Independent',
    ],
    python_requires='>=3.6',
    install_requires=[
        # List your Python dependencies here
    ],
    cmdclass={
        'install': InstallWithEOBLib,
    },
)
