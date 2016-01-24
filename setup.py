from setuptools import setup

def readme():
    with open('README.md') as f:
        return f.read()

setup(name='positronium',
      version='0.1.2',
      description='python tools pertaining to positronium',
      url='https://github.com/PositroniumSpectroscopy/positronium',
      author='Adam Deller',
      author_email='adam.deller1@gmail.com',
      license='MIT',
      packages=['positronium'],
      install_requires=[
          'scipy',
      ],
      include_package_data=True,
      zip_safe=False)