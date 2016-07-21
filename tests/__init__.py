# -*- coding: utf-8 -*-
from biseqt.util import ProgressIndicator
from biseqt.util import Logger

# Turn off all logging for tests
ProgressIndicator.write = lambda *args, **kwargs: None
Logger.log = lambda *args, **kwargs: None
