import setuptools

setuptools.setup(
    name="starpa",
    version="0.0.1",
    url="https://github.com/luidale/starpa",

    author="Hannes Luidalepp",
    author_email="luidale@gmail.com",

    description="Stable RNA processing product analyzer",
    long_description=open('README.rst').read(),

    packages=setuptools.find_packages(),

    install_requires=[],

    classifiers=[
        'Development Status :: 2 - Pre-Alpha',
        'Programming Language :: Python',
        'Programming Language :: Python :: 2',
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.4',
        'Programming Language :: Python :: 3.5',
    ],
)
