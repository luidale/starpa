#import src.starpa

import setuptools
import src.starpa

setuptools.setup(
    name="starpa",
    version=src.starpa.__version__,
    url="https://github.com/luidale/starpa",
	download_url="https://github.com/luidale/starpa/archive/v0.1-alpha.tar.gz",
    author="Hannes Luidalepp",
    author_email="luidale@gmail.com",

    description="Stable RNA processing product analyzer",
    long_description=open('README.rst').read(),

    packages=setuptools.find_packages('src'),
	package_dir={'': 'src'},
	
    install_requires=["pyfaidx","schema","docopt"],

    classifiers=[
        'Development Status :: 2 - Pre-Alpha',
		'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.4',
        'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: 3.6',
		'Environment :: Console',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
        'Operating System :: Unix',
        'Topic :: Scientific/Engineering',
		'Topic :: Scientific/Engineering :: Bio-Informatics'
    ],
    package_data={'': ['data/*.*']}
)
