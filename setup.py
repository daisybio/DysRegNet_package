from setuptools import setup, find_packages
import pathlib


# The directory containing this file
HERE = pathlib.Path(__file__).parent


# The text of the README file
README = (HERE / "README.md").read_text()


setup(name='dysregnet',
      version='0.0.1',
      description='DysRegNet',
      long_description=README,
      long_description_content_type="text/markdown",
      url='',
      author='Zakaria Louadi, olga lazareva ',
      author_email='zakaria.louadi@tum.de, olga.lazareva@tum.de',
      license='GPLv3',
      classifiers=[
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
      ],
      packages=find_packages(),
      include_package_data=True,
      python_requires='>=3.6',
      install_requires=[
        'pandas',
        'numpy>= 1.19',
        'plotly>=4.14.3',
        'scipy>=1.6.2',
        'statsmodels>=0.12.2',
        'tqdm',
         'sklearn', 
    ],

)