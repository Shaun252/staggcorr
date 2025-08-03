from setuptools import setup, find_packages

with open("install_requires.txt") as f:
    dependencies = f.read().strip().split('\n')

setup(name='stagcorr',
      version='0.1.0',
      description='Staggered fermion correlation function solver for lattice QCD',
      url='https://github.com/shaun252/stagcorr',
      author='Shaun Lahert',
      author_email='shaun.lahert@gmail.com',
      license='MIT',
      packages=find_packages(),
      install_requires=dependencies,
      python_requires='>=3.8',
      zip_safe=False)