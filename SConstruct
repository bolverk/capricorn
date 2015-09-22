import glob
import os
from argparse import ArgumentParser

#parser = ArgumentParser(description='gets name of input and output files')
#parser.add_argument('main_source_file',help='path to main source file')
#parser.add_argument('output_file',help='name of output file')
#args = parser.parse_args()

debug = ARGUMENTS.get('debug',0)
source = ARGUMENTS.get('source',None)
target = ARGUMENTS.get('target',None)
compiler = ARGUMENTS.get('compiler','g++')

if compiler=='g++':
    cflags = '-Wfatal-errors'
    if int(debug):
        cflags +=' -O0 -g -pg'
    else:
        cflags +=' -O3'
elif compiler=='clang++':
    cflags = '-Weverything -Werror -ferror-limit=1 -Wno-error=padded'
    if int(debug):
        cflags += ' -O0 -g -pg'
    else:
        cflags += ' -O3'
env = Environment(ENV = os.environ,
                  CXX=compiler,
                  CPPPATH='source',
                  LIBPATH=['.',os.environ['HDF5_LIB_PATH']],
                  LIBS=['capricorn','hdf5','hdf5_cpp'],
                  CXXFLAGS=cflags)
lib_file = env.Library('capricorn',glob.glob('source/*.cpp'))
if None!=source and None!=target:
    env.Program(target,[source])


