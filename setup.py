import src.starpa

import setuptools

setuptools.setup(
    name="starpa",
    version=starpa.__version__,
    url="https://github.com/luidale/starpa",

    author="Hannes Luidalepp",
    author_email="luidale@gmail.com",

    description="Stable RNA processing product analyzer",
    long_description=open('README.rst').read(),

    packages=setuptools.find_packages(),

    install_requires=["pyfaidx","schema"],

    classifiers=[
        'Development Status :: 2 - Pre-Alpha',
		'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.4',
        'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: 3.6',
		'Environment :: Console',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
        'Operating System :: OS Independent',
        'Topic :: Scientific/Engineering',
		'Topic :: Scientific/Engineering :: Bio-Informatics'
    ],
    package_data={'': ['data/*.*','data/tests/*.*']}
)
