#import src.starpa

import setuptools
import src.starpa

setuptools.setup(
    name="starpa",
    version=src.starpa.__version__,
    url="https://github.com/luidale/starpa",
	download_url="https://github.com/luidale/starpa/archive/v0.3.0.tar.gz",
    author="Hannes Luidalepp",
    author_email="luidale@gmail.com",
	license=src.starpa.__license__,
    description="Stable RNA processing product analyzer",
    long_description=open('README.rst').read(),

    packages=setuptools.find_packages('src'),
	package_dir={'': 'src'},
	entry_points = {
        'console_scripts': ['starpa=starpa.__main__:main']},
		
    install_requires=["pyfaidx","schema","docopt","cutadapt"],

    classifiers=[
        'Development Status :: 4 - Beta',
		'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.4',
        'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: 3.6',
		'Programming Language :: Python :: Implementation :: PyPy3',
		'Environment :: Console',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
        'Operating System :: Unix',
        'Topic :: Scientific/Engineering',
		'Topic :: Scientific/Engineering :: Bio-Informatics'
    ],
    package_data={'': ['data/*.*']}
)
