import glob
import os

mode = ARGUMENTS.get('mode','gcc_release')

if mode=='gcc_release':
    compiler = 'g++'
    cflags = ' -Wfatal-errors -O3 '
elif mode=='gcc_debug':
    compiler = 'g++'
    cflags = ' -Wfatal-errors -O0 -g -pg '
elif mode=='clang_release':
    compiler= 'clang++'
    cflags = '-Weverything -Werror -ferror-limit=1 -Wno-error=padded -O3 '
elif mode=='parallel_release':
    compiler='mpiCC'
    cflags = ' -O3 -DWITH_MPI -Wall -Wextra -pedantic -Wfatal-errors '
else:
    raise NameError('unsupported mode')

build_dir = 'build/'+mode
env = Environment(ENV = os.environ,
                  CXX=compiler,
                  CPPPATH='source',
                  LIBPATH=['.',os.environ['HDF5_LIB_PATH']],
                  LIBS=['hdf5','hdf5_cpp'],
                  CXXFLAGS=cflags)
env.VariantDir(build_dir,'source')
lib_file = env.Library(build_dir+'/capricorn',
                       Glob(build_dir+'/*.cpp'))


