# -*- coding: utf-8 -*-
from StringIO import StringIO
from biseqt.util import ProgressIndicator


def test_progress_indicator():
    logs = StringIO()
    ProgressIndicator.write = lambda self, message: logs.write(message)

    indic = ProgressIndicator(num_total=1)
    indic.start()
    indic.progress()
    assert logs.getvalue().strip() == '0/1 \r1/1', \
        'Counting progress indicator works'

    logs = StringIO()
    indic = ProgressIndicator(num_total=1, percentage=True)
    indic.start()
    indic.progress()
    assert logs.getvalue().strip() == '0% \r100%', \
        'Percentage progress indicator works'
