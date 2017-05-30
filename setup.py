from setuptools import setup, find_packages

setup(
    name='mg-process-fastq',
    packages=find_packages(),
    include_package_data=True,
    install_requires=[
        'numpy', 'h5py', 'scipy', 'matplotlib', 'pysam', 'MACS2', 'rpy2',
        'pytest'
    ],
    setup_requires=[
        'pytest-runner',
    ],
    tests_require=[
        'pytest',
    ],
)
