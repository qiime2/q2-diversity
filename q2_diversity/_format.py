# ----------------------------------------------------------------------------
# Copyright (c) 2016-2021, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import qiime2.plugin.model as model

class ProcrustesM2StatisticFmt(model.TextFileFormat):
    def validate(*args):
        pass

ProcrustesM2StatisticDirFmt = model.SingleFileDirectoryFormat(
    'ProcrustesM2StatisticDirFmt', 'ProcrustesM2Statistic.csv', 
    ProcrustesM2StatisticFmt)
