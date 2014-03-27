
def configuration(parent_package='',top_path=None):
   from numpy.distutils.misc_util import Configuration
   config = Configuration('recon3d', parent_package, top_path)
   config.add_extension('_reconutils',sources=['_reconutils.f90'],
            libraries=['lapack','refblas','gfortran'],
            library_dirs=[],
            extra_f90_compile_args=['-O3'])
#   config.add_data_dir('./test_data')
   return config

if __name__ == '__main__':
   from distutils.dir_util import remove_tree
   from numpy.distutils.core import setup
   if os.path.exists('./build'):  
      remove_tree('./build')
   setup(**configuration(top_path='').todict())

