# install specific version via: pip install cffi==0.8.2
# where 0.8.2 comes from: apt-cache policy python-cffi (which is the "backend")
import cffi
import os
import sys
import re

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


# FIXME this should be moved to seq.py
def hp_tokenize(string):
    """Generates (yields) homopolymeric stretches of the given sequences in
    order in tuples of the form ``(char, num, pos)``. For example::

        hp_tokenize('AAACCG') #=> [('A', 3, 0), ('C', 2, 3), ('G', 1, 5)]
    """
    for match in re.finditer(r'(.)\1*', string):
        match, pos = match.group(0), match.start()
        yield match[0], len(match), pos


class ProgressIndicator(object):
    def __init__(self, msg, num_total, percentage=True):
        self.msg, self.num_total = msg, int(num_total)
        self.freq = max(self.num_total/100, 1)
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
