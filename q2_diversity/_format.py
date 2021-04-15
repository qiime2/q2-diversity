# ----------------------------------------------------------------------------
# Copyright (c) 2016-2021, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import qiime2.plugin.model as model
from qiime2.plugin import ValidationError
import qiime2


class ProcrustesM2StatisticFmt(model.TextFileFormat):
    METADATA_COLUMNS = {
        'true M^2 value',
        'p-value for true M^2 value',
        'number of Monte Carlo permutations',
    }

    def validate(self, level):
        try:
            md = qiime2.Metadata.load(str(self))
        except qiime2.metadata.MetadataFileError as md_exc:
            raise ValidationError(md_exc) from md_exc

        for column in sorted(self.METADATA_COLUMNS):
            try:
                md.get_column(column)
            except ValueError as md_exc:
                raise ValidationError(md_exc) from md_exc


ProcrustesM2StatDFmt = model.SingleFileDirectoryFormat(
    'ProcrustesM2StatDFmt', 'ProcrustesM2Statistic.csv',
    ProcrustesM2StatisticFmt)
