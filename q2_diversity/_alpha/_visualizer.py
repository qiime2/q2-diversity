# ----------------------------------------------------------------------------
# Copyright (c) 2016--, QIIME development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import json
import os
import shutil

import scipy
import numpy
import qiime
import pandas as pd
from statsmodels.sandbox.stats.multicomp import multipletests
from trender import TRender


def alpha_group_significance(output_dir: str, alpha_diversity: pd.Series,
                             metadata: qiime.Metadata) -> None:
    metadata_df = metadata.to_dataframe()
    metadata_df = metadata_df.select_dtypes(exclude=[numpy.number])
    categories = metadata_df.columns

    filenames = []
    for category in categories:
        metadata_category = metadata.get_category(category).to_series()
        filtered_metadata = metadata_category[alpha_diversity.index]
        data = pd.concat([alpha_diversity, filtered_metadata], axis=1)
        names = []
        groups = []
        for name, group in data.groupby(filtered_metadata.name):
            names.append('%s (n=%d)' % (name, len(group)))
            groups.append(list(group[alpha_diversity.name]))

        if (len(groups) > 1 and len(groups) != len(data)):
            filename = 'category-%s.jsonp' % category
            filenames.append(filename)
            df = pd.Series(groups, index=names)

            # perform Kruskal-Wallis across all groups
            kw_H_all, kw_p_all = scipy.stats.mstats.kruskalwallis(*groups)
            kwAll = json.dumps({'H': kw_H_all, 'p': kw_p_all})

            # perform pairwise Kruskal-Wallis across all pairs of groups and
            # correct for multiple comparisons
            kw_H_pairwise = []
            for i in range(len(names)):
                for j in range(i):
                    H, p = scipy.stats.mstats.kruskalwallis(groups[i],
                                                            groups[j])
                    kw_H_pairwise.append([names[j], names[i], H, p])
            kw_H_pairwise = pd.DataFrame(
                kw_H_pairwise, columns=['Group 1', 'Group 2', 'H', 'p-value'])
            kw_H_pairwise.set_index(['Group 1', 'Group 2'], inplace=True)
            kw_H_pairwise['q-value'] = multipletests(
                kw_H_pairwise['p-value'], method='fdr_bh')[1]
            kw_H_pairwise.sort_index(inplace=True)
            outfile = 'kruskal-wallis-pairwise-%s.csv' % category
            outpath = os.path.join(output_dir, outfile)
            kw_H_pairwise.to_csv(outpath)

            with open(os.path.join(output_dir, filename), 'w') as fh:
                fh.write('load_data("%s",`' % category)
                df.to_json(fh, orient='split')
                fh.write("`,'%s',`" % kwAll)
                kw_H_pairwise.to_html(fh, classes="table table-striped "
                                                  "table-hover")
                fh.write("`,'./%s');" % outfile)

    TEMPLATES = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                             'assets')
    index = TRender('index.template', path=TEMPLATES)
    rendered_index = index.render({'categories': filenames})
    with open(os.path.join(output_dir, 'index.html'), 'w') as fh:
        fh.write(rendered_index)

    shutil.copytree(os.path.join(TEMPLATES, 'dst'),
                    os.path.join(output_dir, 'dist'))
