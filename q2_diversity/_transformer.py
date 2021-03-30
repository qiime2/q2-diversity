# ----------------------------------------------------------------------------
# Copyright (c) 2016-2021, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import pandas as pd
import qiime2

from .plugin_setup import plugin
from ._format import ProcrustesM2StatisticFmt


@plugin.register_transformer
def _1(data: pd.DataFrame) -> ProcrustesM2StatisticFmt:
    ff = ProcrustesM2StatisticFmt()
    qiime2.Metadata(data).save(str(ff))
    return ff


@plugin.register_transformer
def _2(ff: ProcrustesM2StatisticFmt) -> pd.DataFrame:
    return qiime2.Metadata.load(str(ff)).to_dataframe()


@plugin.register_transformer
def _3(ff: ProcrustesM2StatisticFmt) -> qiime2.Metadata:
    return qiime2.Metadata.load(str(ff))
