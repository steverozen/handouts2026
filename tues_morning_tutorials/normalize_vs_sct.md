# `SCTransform()` vs `NormalizeData` → `FindVariableFeatures` → `ScaleData`

What actually differs in the final Seurat object.

## What ends up in the object

### Standard workflow (`NormalizeData` → `FindVariableFeatures` → `ScaleData`)

- Everything stays in the `RNA` assay. `DefaultAssay` is unchanged.
- `layers$counts` untouched, `layers$data` = `log1p(counts / colSums * 10000)`, `layers$scale.data` = z-scored `data`, restricted to the 2000 variable features and clipped at ±10.
- `VariableFeatures()` = 2000 genes ranked by VST-standardized variance.
- `meta.features` (or `[[assay]][[]]` in v5) gains `vf_vst_counts_mean`, `vf_vst_counts_variance`, `vf_vst_counts_variance.standardized`.
- No new `meta.data` columns.

### `SCTransform()`

- Creates a **new assay** `SCT` and sets `DefaultAssay(obj) <- "SCT"`. The `RNA` assay stays, with counts only and no `data` layer, so `RNA` is left unnormalized unless you also ran `NormalizeData`.
- `SCT@counts` = *corrected* counts (counts reversed out of the model as if every cell had the same sequencing depth), `SCT@data` = `log1p` of those corrected counts, `SCT@scale.data` = **Pearson residuals**, not z-scores, clipped to ±`sqrt(ncol/30)` by default and centered but not variance-scaled.
- `SCT@scale.data` holds only the variable features by default (`return.only.var.genes = TRUE`), so residuals for other genes are not stored.
- `VariableFeatures()` = **3000** genes by default, ranked by *residual variance*, a different criterion and a different gene list than VST.
- **Fewer genes**: features detected in fewer than `min_cells = 5` cells are dropped, so `rownames(obj[["SCT"]])` is a strict subset of `rownames(obj[["RNA"]])`.
- A fitted model is stored: `obj[["SCT"]]@SCTModel.list`, with per-gene `theta`, `(Intercept)`, `log_umi`, gmean, detection rate, residual variance, plus the cell attributes and `umi.assay` / `min_variance` / `clip.range` settings. Nothing comparable exists in the standard path.
- `meta.data` gains `nCount_SCT` and `nFeature_SCT`.

## Side by side

| | Standard | `SCTransform` |
|---|---|---|
| Assay written | `RNA` | new `SCT`, becomes default |
| `data` layer | `log1p` CP10K | `log1p` of corrected counts |
| `scale.data` layer | z-scores, clipped ±10 | Pearson residuals, clipped ±`sqrt(n/30)` |
| Variable features | 2000, VST variance | 3000, residual variance |
| Gene set | all genes | drops genes in < 5 cells |
| Depth correction | one global size factor | per-gene NB GLM on `log_umi` |
| Model stored | none | `SCTModel.list` |
| New `meta.data` | none | `nCount_SCT`, `nFeature_SCT` |

## Practical consequences

- **Depth correction**: SCT regresses sequencing depth in a per-gene NB GLM. The log-normalize path applies one global size factor, which under-corrects highly expressed genes and over-corrects lowly expressed ones.
- **`vars.to.regress`** goes into `SCTransform()` rather than `ScaleData()`. Regressing `percent.mt` there is the SCT equivalent.
- **Downstream DE**: `FindMarkers()` on the `SCT` assay uses the `data` layer. If you ran `SCTransform` separately per sample and then merged or integrated, you must call `PrepSCTFindMarkers()` first, otherwise the models are not on a common scale. No such step exists for log-normalized data.
- **Visualization**: `FeaturePlot` / `VlnPlot` will silently switch to `SCT` values after `SCTransform` because the default assay changed. Many people prefer to plot from `RNA` with `NormalizeData` also run, which is why the common recipe is to run both.
- **Integration**: use `normalization.method = "SCT"` in `FindIntegrationAnchors` / `IntegrateData`, or `SCTransform` + `IntegrateLayers(..., normalization.method = "SCT")` in v5.
- **`vst.flavor = "v2"`** is the default in recent versions, which changes the theta regularization relative to older SCT results, so objects built by different Seurat versions are not directly comparable.

## References

- Hafemeister C, Satija R (2019). Normalization and variance stabilization of single-cell RNA-seq data using regularized negative binomial regression. *Genome Biology* 20:296. https://doi.org/10.1186/s13059-019-1874-1
- Choudhary S, Satija R (2022). Comparison and evaluation of statistical error models for scRNA-seq. *Genome Biology* 23:27. https://doi.org/10.1186/s13059-021-02584-9
