import qiime2.plugin.model as model

class M2CalcFmt(model.TextFileFormat):
    def validate(*args):
        pass

M2CalcDirFmt = model.SingleFileDirectoryFormat(
    'M2CalcDirFmt', 'm2.csv', M2CalcFmt)
