from setuptools import setup, find_packages

setup(
    name='mg-process-fastq',
    packages=find_packages(),
    include_package_data=True,
    install_requires=[
<<<<<<< HEAD
        'numpy', 'h5py', 'scipy', 'matplotlib', 'pysam==0.9.1.4', 'MACS2',
        'rpy2', 'pytest'
=======
#        'numpy', 'h5py', 'scipy', 'matplotlib', 'pysam', 'MACS2', 'rpy2',
#        'pytest'
>>>>>>> rf-tests-MNase
    ],
    setup_requires=[
        'pytest-runner',
    ],
    tests_require=[
        'pytest',
    ],
)
