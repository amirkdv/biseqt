#!/usr/env/bin python

# install specific version via: pip install cffi==0.8.2
# where 0.8.2 comes from: apt-cache policy python-cffi (which is the "backend")
import cffi
import os

ffi = cffi.FFI()
lib = ffi.dlopen(os.path.join(os.path.dirname(__file__), 'libalign.so'))
with open(os.path.join(os.path.dirname(__file__), 'libalign.h')) as f:
    ffi.cdef(f.read())
