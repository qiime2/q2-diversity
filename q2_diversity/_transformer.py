# format from q2-dada2 and q2-deblur 
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
