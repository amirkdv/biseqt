# install specific version via: pip install cffi==0.8.2
# where 0.8.2 comes from: apt-cache policy python-cffi (which is the "backend")
import cffi
import os

ffi = cffi.FFI()
lib = ffi.dlopen(os.path.join(os.path.dirname(__file__), 'libalign.so'))
with open(os.path.join(os.path.dirname(__file__), 'libalign.h')) as f:
    ffi.cdef(f.read())

class CffiObject(object):
    """Generic cffi wrapper for C structs, delegates all unknown attributes to
    the underlying C pointer. Subclasses must populate self.c_obj in their
    constructors with a pointer to their underlying C struct.

    Attributes:
        c_obj (cffi.cdata): the underlying C pointer.
    """
    def __init__(self, c_type, **kw):
        self.c_obj = ffi.new('%s *' % c_type)

    def __getattr__(self, name):
        return getattr(self.c_obj, name)

