from setuptools import setup, find_packages

packages = find_packages(".")

setup(
    name="acts",
    version="0.1",
    description="ACTS",
    packages=packages,
    install_requires=[],
    package_data = {
        'acts': ["*.so"]
    },
    include_package_data=True
)