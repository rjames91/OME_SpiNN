try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup

from collections import defaultdict
import os


# Build a list of all project modules, as well as supplementary files
main_package = "spinnak_ear"
extensions = {".aplx", ".boot", ".cfg", ".json", ".sql", ".template", ".xml",
              ".xsd"}
main_package_dir = os.path.join(os.path.dirname(__file__), main_package)
start = len(main_package_dir)
packages = []
package_data = defaultdict(list)
for dirname, dirnames, filenames in os.walk(main_package_dir):
    if '__init__.py' in filenames:
        package = "{}{}".format(
            main_package, dirname[start:].replace(os.sep, '.'))
        packages.append(package)
    for filename in filenames:
        _, ext = os.path.splitext(filename)
        if ext in extensions:
            package = "{}{}".format(
                main_package, dirname[start:].replace(os.sep, '.'))
            package_data[package].append(filename)

setup(
    name="SpiNNak-Ear",
    version="1.0",
    description="Ear",
    url="https://github.com/rjames91/OME_SpiNN",
    license="GNU GPLv3.0",
    packages=packages,
    package_data=package_data,
    install_requires=[
        'SpiNNUtilities >= 1!4.0.1, < 1!5.0.0',
        'SpiNNStorageHandlers >= 1!4.0.1, < 1!5.0.0',
        'SpiNNMachine >= 1!4.0.1, < 1!5.0.0',
        'SpiNNMan >= 1!4.0.1, < 1!5.0.0',
        'SpiNNaker_PACMAN >= 1!4.0.1, < 1!5.0.0',
        'SpiNNaker_DataSpecification >= 1!4.0.1, < 1!5.0.0',
        'spalloc >= 1.0.1, < 2.0.0',
        'SpiNNFrontEndCommon >= 1!4.0.1, < 1!5.0.0',
        'numpy', 'lxml', 'six', 'csa']
)
