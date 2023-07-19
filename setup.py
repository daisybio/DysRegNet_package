from setuptools import setup, find_packages
import pathlib


# The directory containing this file
HERE = pathlib.Path(__file__).parent


# The text of the README file
README = (HERE / "README.md").read_text()


setup(name='dysregnet',
      version='0.0.4',
      description='DysRegNet',
      long_description=README,
      long_description_content_type="text/markdown",
      url='https://github.com/biomedbigdata/DysRegNet_package',
      author='Zakaria Louadi, Olga Lazareva, Johannes Kersting',
      author_email='zakaria.louadi@tum.de, olga.lazareva@tum.de, johannes.kersting@tum.de',
      license='GPLv3',
      classifiers=[
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
      ],
      package_dir = {'': 'src'},
      packages=['dysregnet'],
      include_package_data=True,
      python_requires='>=3.7',
      install_requires=[
        'pandas',
        'numpy>= 1.19',
        'scipy',
        'statsmodels',
        'tqdm',
        'scikit-learn', 

    ],

)
