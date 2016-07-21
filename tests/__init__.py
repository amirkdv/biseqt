# -*- coding: utf-8 -*-
from biseqt import ProgressIndicator
from biseqt import Logger

# skip all logging
ProgressIndicator.write = lambda *args: None
Logger.log = lambda *args, **kwargs: None
