---
name: Core diversity metrics.
description: Applies a collection of diversity metrics to a feature table.
inputs:
    - table:
        - FeatureTable[Frequency]
        - biom.Table
    - phylogeny:
        - Phylogeny[Rooted]
        - skbio.TreeNode
parameters:
    - sampling_depth:
        - Int
        - int
outputs:
    - faith_pd_vector:
        - SampleData[AlphaDiversity]
        - pandas.Series
    - observed_otus_vector:
        - SampleData[AlphaDiversity]
        - pandas.Series
    - shannon_vector:
        - SampleData[AlphaDiversity]
        - pandas.Series
    - evenness_vector:
        - SampleData[AlphaDiversity]
        - pandas.Series
    - unweighted_unifrac_distance_matrix:
        - DistanceMatrix
        - skbio.DistanceMatrix
    - weighted_unifrac_distance_matrix:
        - DistanceMatrix
        - skbio.DistanceMatrix
    - jaccard_distance_matrix:
        - DistanceMatrix
        - skbio.DistanceMatrix
    - bray_curtis_distance_matrix:
        - DistanceMatrix
        - skbio.DistanceMatrix
    - unweighted_unifrac_pcoa_results:
        - PCoAResults
        - skbio.OrdinationResults
    - weighted_unifrac_pcoa_results:
        - PCoAResults
        - skbio.OrdinationResults
    - jaccard_pcoa_results:
        - PCoAResults
        - skbio.OrdinationResults
    - bray_curtis_pcoa_results:
        - PCoAResults
        - skbio.OrdinationResults
---
## Core diversity metrics

This method rarefies a feature table to an even sampling depth and then applies a suite of alpha and beta diversity metrics. The metrics applied are:
 * Alpha diversity
  * Faith's Phylogenetic Diversity index: a qualitative phylogenetic measure of community richness
  * Observed OTUs: a qualitative non-phylogenetic measure of community richness
  * Shannon's Diversity index: a quantitative non-phylogenetic measure of community richness
  * Evenness (or Pielou's Evenness): a non-phylogenetic measure of community evenness
 * Beta diversity
  * Unweighted UniFrac: a qualitative phylogenetic measure of pairwise community dissimilarity (computed for all pairs of samples)
  * Weighted UniFrac: a quantitative phylogenetic measure of pairwise community dissimilarity (computed for all pairs of samples)
  * Jaccard: a qualitative non-phylogenetic measure of pairwise community dissimilarity (computed for all pairs of samples)
  * Bray-Curtis: a quantitative non-phylogenetic measure of pairwise community dissimilarity (computed for all pairs of samples)

Principle coordinates analysis is subsequently applied to all pairwise distance matrices resulting from the beta diversity computations.

### Rarefy feature table

```python
>>> from q2_feature_table import rarefy
>>> rarefied_table = rarefy(table=table, sampling_depth=sampling_depth)
```

### Compute alpha diversity metrics

```python
>>> from q2_diversity import alpha, alpha_phylogenetic
>>> faith_pd_vector = alpha_phylogenetic(table=rarefied_table, phylogeny=phylogeny, metric='faith_pd')
>>> observed_otus_vector = alpha(table=rarefied_table, metric='observed_otus')
>>> shannon_vector = alpha(table=rarefied_table, metric='shannon')
>>> evenness_vector = alpha(table=rarefied_table, metric='pielou_e')
```

### Compute beta diversity metrics

```python
>>> from q2_diversity import beta, beta_phylogenetic
>>> unweighted_unifrac_distance_matrix = beta_phylogenetic(table=rarefied_table, phylogeny=phylogeny, metric='unweighted_unifrac')
>>> weighted_unifrac_distance_matrix = beta_phylogenetic(table=rarefied_table, phylogeny=phylogeny, metric='weighted_unifrac')
>>> jaccard_distance_matrix = beta(table=rarefied_table, metric='jaccard')
>>> bray_curtis_distance_matrix = beta(table=rarefied_table, metric='braycurtis')
```

### Apply principal coordinate analysis to beta diversity results

```python
>>> from q2_diversity import pcoa
>>> unweighted_unifrac_pcoa_results = pcoa(distance_matrix=unweighted_unifrac_distance_matrix)
>>> weighted_unifrac_pcoa_results = pcoa(distance_matrix=weighted_unifrac_distance_matrix)
>>> jaccard_pcoa_results = pcoa(distance_matrix=jaccard_distance_matrix)
>>> bray_curtis_pcoa_results = pcoa(distance_matrix=bray_curtis_distance_matrix)
```
