import cffi
import os
import sys

ffi = cffi.FFI()
lib = ffi.dlopen(os.path.join(os.path.dirname(__file__), 'pwlib', 'pwlib.so'))
with open(os.path.join(os.path.dirname(__file__), 'pwlib', 'pwlib.h')) as f:
    ffi.cdef(f.read())
    # so we can call strlen from python:
    # FIXME used only once in pw.py, can we get rid of it?
    ffi.cdef('size_t strlen(const char*);')


class CffiObject(object):
    """Generic cffi wrapper for C structs, delegates all unknown attributes to
    the underlying C pointer. Subclasses must populate :attr:`c_obj` in their
    constructors with a pointer to their underlying C struct.

    Attributes:
        c_obj (cffi.cdata): the underlying C pointer.
    """
    def __init__(self, c_type, **kw):
        self.c_obj = ffi.new('%s *' % c_type)

    def __getattr__(self, name):
        return getattr(self.c_obj, name)


class ProgressIndicator(object):
    def __init__(self, msg, num, percentage=True):
        self.msg, self.num = msg, int(num)
        self.freq = max(self.num/100, 1)
        self.progress_cnt = 0
        self.percentage = percentage
        if self.percentage and self.num == 0:
            raise ValueError('Percentage progress indicator needs > 0 rounds')

    def start(self):
        self.status()

    def finish(self):
        self.status()
        sys.stderr.write('.\n')

    def status(self):
        if self.percentage:
            sys.stderr.write('\r%s: %d%%' %
                             (self.msg, (self.progress_cnt * 100) / self.num))
        else:
            sys.stderr.write('\r%s: %d/%d' %
                             (self.msg, self.progress_cnt, self.num))

    def progress(self, num=1):
        self.progress_cnt += num
        if not self.percentage or self.progress_cnt % self.freq == 0:
            self.status()
