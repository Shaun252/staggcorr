from setuptools import setup, find_packages
from itertools import product

data_dirs = ['stagcorr']

include_subdirs = ['correlators',
                   'gauge'
                   'lattice'
                   'propagatorUtils']

suffixes = ['.py']

data_paths = ["{}/{}/*{}".format(dr, subdir, suffix)
              for dr, subdir, suffix in product(data_dirs, include_subdirs,
                                                suffixes)]


print(find_packages())

with open("install_requires.txt") as f:
    dependencies = f.read().split()

setup(name='stagcorr',
      version='0.1',
      description='staggered correlation function solver',
      url='http://github.com/shaun252',
      author='Shaun Lahert',
      author_email='slahert@illinois.edu',
      license='NA',
      packages=find_packages(),
      package_dir={'': '.'},
      package_data={'stagcorr': data_paths},
      install_requires=dependencies,
      zip_safe=False)