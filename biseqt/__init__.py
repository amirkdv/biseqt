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
    def __init__(self, quiet=False, f=sys.stderr, num_total=None, percentage=False):
        self.f = f
        self.quiet = quiet
        self.percentage = percentage
        if self.percentage:
            assert isinstance(num_total, int) and num_total > 0
            self.num_total = num_total
            self.freq = max(self.num_total/100, 1)
        else:
            self.num_total = num_total
            self.freq = None

        self.progress_cnt = 0

    def start(self):
        if self.quiet:
            self
        self.status()

    def finish(self):
        if self.quiet:
            return
        self.status()
        self.f.write('\r')

    def status(self):
        if self.quiet:
            return
        if self.percentage:
            self.f.write('\r%d%% ' % (self.progress_cnt * 100. / self.num_total))
        else:
            num = str(self.num_total) if self.num_total else ''
            self.f.write('\r%d/%s ' % (self.progress_cnt, num))

    def progress(self, num=1):
        self.progress_cnt += num
        if not self.percentage or self.progress_cnt % self.freq == 0:
            self.status()
