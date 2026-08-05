# Pseudobulk DESeq2 when a cluster has very few cells per sample

**Situation:** pseudobulk differential expression for one cluster (dendritic
cells) across affected versus control samples, where the number of dendritic
cells contributing to each sample's pseudobulk profile ranges from 2 to 34.

## Short answer

No, do not put all samples into the DESeq2 fit as they stand. A pseudobulk
profile built from 2 cells is not a noisy measurement of that sample's
dendritic cell state, it is a measurement of 2 particular cells. Apply a
pre-specified minimum cell count per sample (10 cells is the usual default,
20 is defensible), and a minimum library size, drop the samples that fail,
and then check whether enough samples survive in each group to run the test
at all. With a maximum of 34 cells anywhere in the dataset, expect that this
cluster may simply not support a reliable DE analysis, and be prepared to
report that rather than to report a gene list.

## Why very low cell counts break the pseudobulk assumptions

Pseudobulk works because summing counts over hundreds of cells averages away
cell-level sampling noise and dropout, leaving a per-sample profile whose
residual variability is dominated by biological differences between
individuals. That is the variability DESeq2's negative binomial dispersion is
meant to capture.

With 2 to 34 cells that averaging does not happen.

1. **Cell-level sampling dominates.** The pseudobulk profile from 2 cells is
   driven by which 2 cells were captured, including their cell cycle phase,
   ambient RNA contamination, and any doublet or misclassification. One
   mislabeled cell can move the whole sample.

2. **Heteroscedasticity that the model cannot represent.** DESeq2 fits one
   dispersion per gene, shared across all samples, after size factor
   normalization. A 2 cell sample and a 34 cell sample have wildly different
   true precision. The shared dispersion gets inflated to accommodate the
   worst samples, which costs power for every gene, including in the samples
   that were fine.

3. **Sparsity and size factors.** Tiny pseudobulk libraries have many zero
   and near zero genes. Median of ratios size factor estimation degrades when
   a sample shares few nonzero genes with the rest, and normalization errors
   then propagate into the log fold changes.

4. **Outlier driven calls.** Cook's distance filtering in DESeq2 helps, but
   with small group sizes it either fails to flag the offending sample or
   removes so many genes that the analysis becomes unstable.

## The threshold to use, and how to choose it

Standard practice in the pseudobulk literature is to require a minimum number
of cells per sample per cluster before that sample contributes. The muscat
workflow (Crowell et al., 2020) filters on both cell number and expression,
and 10 cells per sample is the common default. Some groups use 20.

Two rules matter more than the exact number:

- **Fix the threshold before looking at the results.** Choosing 10 versus 20
  after seeing which one gives a nicer volcano plot is a selection procedure
  that invalidates the p values.
- **Filter on total counts as well as on cell number.** Cells differ enormously
  in UMI depth. A sample with 12 low depth dendritic cells can carry fewer
  counts than a sample with 5 deep ones. Require, say, at least 10 cells
  *and* a pseudobulk library size not far below the rest of the retained
  samples.

## What to do before deciding

Build this table first and look at it, one row per sample:

| sample | group | n cells in cluster | pseudobulk total UMI | n genes with count > 0 |

Then:

1. **Check whether cell number is confounded with group.** If the affected
   samples systematically have fewer dendritic cells than controls, that
   difference in abundance is itself a result, and it should be reported
   separately as a compositional finding. It also means that filtering on
   cell number preferentially removes one group, which biases the remaining
   comparison. Say so explicitly in the writeup.

2. **Count survivors per group.** After filtering you want at least 3 samples
   per group, and realistically 4 or 5 for anything you intend to believe.
   Below 3 per group, DESeq2 will still return a table, and that table should
   not be interpreted.

3. **Look at a PCA or sample-sample correlation heatmap of the retained
   pseudobulk profiles.** If the retained samples separate by cell number
   rather than by group, the filter was not strict enough.

## If filtering leaves too few samples

Ranked from most to least preferable.

- **Report that this cluster is underpowered and stop.** This is a legitimate
  and common outcome. Report the per sample cell counts and the abundance
  comparison, and do not report a gene list.

- **Merge the cluster upward.** If the dendritic cells split into subtypes,
  or sit next to a related myeloid population, analyze at a coarser level of
  annotation where each sample contributes enough cells. State that the
  analysis is at the coarser level.

- **Downweight instead of dropping.** `limma-voom` with
  `voomWithQualityWeights()` estimates a weight per sample and lets noisy
  low cell samples contribute proportionally less rather than not at all.
  This keeps n up, which matters when n is the binding constraint. It is a
  reasonable alternative to a hard threshold, and it uses the same pseudobulk
  count matrix, so you can run it alongside DESeq2 as a sensitivity check.

- **Model cells directly with a per sample random effect.** A mixed model on
  cell level counts (for example NEBULA, or `glmmTMB` with a random intercept
  for sample) uses the individual cells while still treating the sample as
  the unit of replication, which is what avoids pseudoreplication bias
  (Zimmerman et al., 2021). This is more work and more fragile with 2 cells
  in a sample, and it does not create information that is not there, but it
  does not throw samples away.

What **not** to do: run cell level Wilcoxon or the default `Seurat::FindMarkers`
across affected versus control cells pooled over samples. That treats cells as
independent replicates and produces large numbers of false positives, which is
the central result of Squair et al., 2021.

## Suggested concrete procedure

1. Aggregate raw counts per sample for the dendritic cell cluster with
   `Seurat::AggregateExpression()` or `muscat::aggregateData()`. Sum raw
   counts, do not average normalized values.
2. Record `n_cells` and pseudobulk library size per sample.
3. Drop samples with `n_cells < 10`, and drop samples whose library size is a
   clear outlier low.
4. If fewer than 3 samples remain in either group, stop and report the cluster
   as underpowered.
5. Gene filter on the retained samples, for example
   `edgeR::filterByExpr()`, or require count >= 10 in at least the size of the
   smaller group.
6. Run DESeq2 with `~ group` plus any needed covariates, and use
   `lfcShrink(type = "apeglm")` for ranking.
7. Report in the methods the threshold used, how many samples were dropped
   from each group, and the per sample cell counts. Put the full table in a
   supplement.
8. Sensitivity check: rerun at a threshold of 20 cells, and rerun with
   `voomWithQualityWeights()` on all samples. If the top genes are unstable
   across those three analyses, the result is not solid.

## References

- Crowell HL, Soneson C, Germain P-L, Calini D, Collin L, Raposo C, Malhotra D,
  Robinson MD. "muscat detects subpopulation-specific state transitions from
  multi-sample multi-condition single-cell transcriptomics data."
  *Nature Communications*, 2020.
  DOI: 10.1038/s41467-020-19894-4
  https://doi.org/10.1038/s41467-020-19894-4

- Squair JW, Gautier M, Kathe C, Anderson MA, James ND, Hutson TH, Hudelle R,
  Qaiser T, Matson KJE, Barraud Q, Levine AJ, La Manno G, Skinnider MA,
  Courtine G. "Confronting false discoveries in single-cell differential
  expression." *Nature Communications*, 2021.
  DOI: 10.1038/s41467-021-25960-2
  https://doi.org/10.1038/s41467-021-25960-2

- Zimmerman KD, Espeland MA, Langefeld CD. "A practical solution to
  pseudoreplication bias in single-cell studies."
  *Nature Communications*, 2021.
  DOI: 10.1038/s41467-021-21038-1
  https://doi.org/10.1038/s41467-021-21038-1
