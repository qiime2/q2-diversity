# format from q2-dada2 and q2-deblur 
import pandas as pd

from .plugin_setup import plugin
from ._format import M2CalcFmt


@plugin.register_transformer
def _1(data: pd.DataFrame) -> M2CalcFmt:
    ff = M2CalcFmt()
    data.to_csv(str(ff))
    return ff


@plugin.register_transformer
def _2(ff: M2CalcFmt) -> pd.DataFrame:
    return pd.read_csv(str(ff), index_col=0)
