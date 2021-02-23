import qiime2.plugin.model as model

class ProcrustesM2StatisticFmt(model.TextFileFormat):
    def validate(*args):
        pass

ProcrustesM2StatisticDirFmt = model.SingleFileDirectoryFormat(
    'ProcrustesM2StatisticDirFmt', 'ProcrustesM2Statistic.csv', 
    ProcrustesM2StatisticFmt)
