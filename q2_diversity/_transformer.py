# ----------------------------------------------------------------------------
# Copyright (c) 2016-2020, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import pandas as pd

from .plugin_setup import plugin
from ._format import ProcrustesM2StatisticFmt


@plugin.register_transformer
def _1(data: pd.DataFrame) -> ProcrustesM2StatisticFmt:
    ff = ProcrustesM2StatisticFmt()
    data.to_csv(str(ff))
    return ff


@plugin.register_transformer
def _2(ff: ProcrustesM2StatisticFmt) -> pd.DataFrame:
    return pd.read_csv(str(ff), index_col=0)
