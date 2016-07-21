# -*- coding: utf-8 -*-
import sys
import logging


class Logger(object):
    logging.basicConfig(
        format='%(levelname)s [%(asctime)s] %(header)s %(message)s'
    )

    def __init__(self, log_level=logging.INFO, header='', f=None):
        """Creates an instance of :class:`Logger` with the given verbosity
        level and optional header."""
        self._logger = logging.getLogger('biseqt')
        self._logger.setLevel(log_level)
        self._header = header

    def log(self, message, level=logging.INFO):
        """Logs a message of given severity level."""
        self._logger.log(level, message, extra={'header': self._header})


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
