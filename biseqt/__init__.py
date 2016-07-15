import cffi
import os
import sys
import logging


class Logger(object):
    logging.basicConfig(
        format='%(levelname)s [%(asctime)s] %(header)s %(message)s'
    )

    def __init__(self, log_level=logging.INFO, header=''):
        """Creates an instance of :class:`Logger` with the given verbosity
        level and optional header."""
        self._logger = logging.getLogger('biseqt')
        self._logger.setLevel(log_level)
        self._header = header

    def log(self, message, level=logging.INFO):
        """Logs a message of given severity level."""
        self._logger.log(level, message, extra={'header': self._header})

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
    """Reports progress to the user:

        >>> import time
        >>> indic = ProgressIndicator(num_total=10)
        >>> indic.start()
        >>> for i in range(10):
        ...     indic.progress()
        ...     time.sleep(0.1)
        >>> indic.finish()

    Attributes:
        f (file): An open file to which updates are writen.
        num_total (int): Total number of expected steps; default is None in
            which case percentage reporting is not allowed and non-percentage
            reporting only shows the current count.
        percentage (bool): Whether to show percentage progress instead of a
            running count; default is False.
    """
    def __init__(self, f=sys.stderr, num_total=None, percentage=False):
        self.f = f
        self.percentage = percentage
        if self.percentage:
            assert isinstance(num_total, int) and num_total > 0
            self.num_total = num_total
            self.freq = max(self.num_total/100, 1)
        else:
            self.num_total = num_total
            self.freq = None

        self.progress_cnt = 0

    def write(self, contents):
        self.f.write(contents)

    def start(self):
        self.status()

    def finish(self):
        self.status()
        self.write('\r')

    def status(self):
        if self.percentage:
            self.write('\r%d%% ' % (self.progress_cnt * 100. / self.num_total))
        else:
            out_of = '/%d' % self.num_total if self.num_total else ''
            self.write('\r%d%s ' % (self.progress_cnt, out_of))

    def progress(self, num=1):
        self.progress_cnt += num
        if not self.percentage or self.progress_cnt % self.freq == 0:
            self.status()
