from distutils.core import setup

setup(name='omsi',
      version='0.1',
      description='OpenMSI Toolkit',
      author='Oliver Ruebel and Ben Bowen',
      author_email=' ',
      url='http://openmsi.nersc.gov',
      packages=['omsi' , 'omsi/analysis' , 'omsi/dataformat' , 'omsi/tools' , 'omsi/viewer'],
      requires=['numpy' , 'h5py']
     )