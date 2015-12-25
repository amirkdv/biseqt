# install specific version via: pip install cffi==0.8.2
# where 0.8.2 comes from: apt-cache policy python-cffi (which is the "backend")
import cffi
import os
import sys

ffi = cffi.FFI()
lib = ffi.dlopen(os.path.join(os.path.dirname(__file__), 'libalign.so'))
with open(os.path.join(os.path.dirname(__file__), 'libalign.h')) as f:
    ffi.cdef(f.read())

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


def hp_tokenize(string):
    """Generator for homopolymeric substrings in a given sequences. Each value
    is a (char, num) tuple, e.g. ``("A", 3)`` means ``AAA``.
    """
    counter = 0
    assert(len(string))
    cur_char = string[0]
    pos = 0 # the starting position of the current h.p. stretch
    for idx, char in enumerate(string):
        if char == cur_char:
            continue
        else:
            # char is the first letter of the next h.p. stretch
            yield cur_char, idx-pos, pos
            pos = idx
            cur_char = char
    # left overs:
    if pos <= len(string):
        yield cur_char, len(string) - pos, pos


class ProgressIndicator(object):
    def __init__(self, msg, num_total, percentage=True):
        self.msg, self.num_total = msg, num_total
        self.freq = max(num_total/100, 1)
        self.progress_cnt = 0
        self.percentage = percentage
        if self.percentage and self.num_total == 0:
            raise ValueError('Cannot make a percentage progress indicator with num_total=0')

    def start(self):
        self.status()

    def finish(self):
        self.status();
        sys.stderr.write('.\n');

    def status(self):
        if self.percentage:
            sys.stderr.write('\r%s: %d%%' %
                (self.msg, (self.progress_cnt * 100)/self.num_total))
        else:
            sys.stderr.write('\r%s: %d/%d' %
                (self.msg, self.progress_cnt, self.num_total))

    def progress(self, num=1):
        self.progress_cnt += num
        if not self.percentage or self.progress_cnt % self.freq == 0:
            self.status()
