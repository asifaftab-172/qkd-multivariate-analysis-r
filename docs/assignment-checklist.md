# M5 Assignment 5.1 requirements and submission checklist

This map uses the assignment wording supplied on September 22, 2026. It records where each requested analysis is addressed and what to verify before submission. Coverage cannot guarantee a particular grade.

## Resolve the task/rubric mismatch through coverage

Task 4 requests **linear discriminant analysis** using `MASS::lda()`. The rubric instead assigns 20 points to **multivariate regression analysis and interpretation**. The project includes both: LDA for classifying five predefined conditions and a multivariate multiple regression for jointly explaining QBER and log secure key rate. Neither section replaces the other. No instructor clarification is assumed.

## Grading criteria: 100 points

| Criterion | Points | Evidence in the project |
| --- | ---: | --- |
| Data preprocessing and visualization | 20 | Source/provenance and variable descriptions; schema, range, and missing-value checks; descriptive summaries; log transformation and standardization; scatter and correlation figures. |
| PCA analysis and interpretation | 20 | Standardized PCA, per-component and cumulative explained variance, stated retention criterion, loadings, scree plot, biplot, and interpretation of the reduced space. |
| Clustering analysis and interpretation | 20 | K-means and Ward hierarchical clustering, elbow and silhouette evidence, dendrogram, cluster PCA plots, quantitative profiles, and comparison with experimental conditions. |
| Multivariate regression analysis and interpretation | 20 | Two-outcome model with attenuation, condition, and their interaction; multivariate tests; outcome-specific coefficients and fit; fitted plots; residual diagnostics and dependence limitations. |
| Overall report clarity and insights | 20 | Structured research question, dataset description, results and interpretation for each method, legible figures/tables, interdisciplinary relevance, ethical and methodological limitations, citations, and a concluding synthesis. |

## Task-by-task coverage

| Task | Required work and evidence |
| --- | --- |
| 1.1 Load the data in R | The submitted R script loads the supplied combined CSV or the five original CSVs. This is an external, real-world experimental dataset; it is not the excluded wine dataset. |
| 1.2 Descriptive EDA | Overall and condition-specific descriptive statistics are written under `results/` and interpreted in the report. |
| 1.3 Missing values | Required fields are checked for missing/non-finite values before analysis, with an explicit explanation of whether imputation is needed. Invalid input produces an error rather than silent deletion. |
| 1.4 Scale as needed | Secure key rate is transformed with `log10(SKR + 1)`; PCA/clustering use standardized quantitative variables. LDA scaling is learned from its training partition. |
| 1.5 Visualize the data | Scatterplot matrix, correlation matrix, and outcome-versus-attenuation figures (`01`–`04`). |
| 2.1 Conduct PCA | PCA of attenuation, QBER, and log secure key rate. |
| 2.2 Choose components | Explicit comparison of retention criteria and the interpretation of the retained component count. |
| 2.3 Explain component variance | Per-component and cumulative variance table and scree plot. |
| 2.4 Biplot and scree plot | `05_pca_scree.pdf` and `06_pca_biplot.pdf`. |
| 2.5 Interpret reduction | Explain the axes using loadings, the information preserved/lost, and the distinction between a two-dimensional display and the retention decision. |
| 3.1 Both clustering methods | K-means and Ward hierarchical clustering of standardized quantitative data. |
| 3.2 Choose cluster count | Elbow/silhouette comparison and hierarchical dendrogram; explain the selected count and any competing solutions. |
| 3.3 Visualize clusters | K-means and hierarchical memberships projected into the PCA space (`08` and `10`), plus dendrogram (`09`). |
| 3.4 Interpret structure | Cluster profiles (`13`) and condition cross-tabulations identify performance regimes and overlap. Unsupervised clusters are interpreted rather than assumed to reproduce the five known classes. |
| 4.1 Define groups/predictors | `NoiseCondition` has Baseline, 3, 7, 9, and 12 dBm groups; quantitative predictors are attenuation, QBER, and log secure key rate. The class-defining noise field is excluded from LDA predictors. |
| 4.2 Fit MASS LDA | `MASS::lda()` is used with a reproducible stratified train/test split. |
| 4.3 Interpret discriminants | Coefficient/loadings table and interpretation of variables separating conditions; signs and correlated predictors are considered. |
| 4.4 Scores and separation plot | Discriminant scores and `11_lda_separation.pdf`; interpret separation and overlap. |
| 4.5 Classification evaluation | Held-out confusion matrix (`12_lda_confusion_matrix.pdf`) with counts and row percentages, overall accuracy, and class-specific metrics. Axis labels identify actual and predicted classes. |
| 4.6 Interpret misclassification | Identify difficult pairs/classes from the actual matrix and qualify the within-campaign evaluation. |
| Rubric regression requirement | `cbind(QBER, LogSKR) ~ VOA_c * NoiseCondition`, multivariate tests, coefficients, fit, `14_regression_fits.pdf`, and `15_regression_diagnostics.pdf`. |
| Interdisciplinary purpose | Connect network engineering with statistical learning, reliability, and secure-communication decisions; discuss implications and the limits of claims. |

The report must interpret the numerical results and diagnostics; merely generating the listed files is not sufficient coverage.

## Required submission files

The automated build places these **three** files in `submission/` and in the downloadable `qkd-analysis-and-submission` workflow artifact:

1. `M5Assign 5.1 Data__Ronggon_Asif.R` — all analysis code.
2. `M5Assign 5.1 Data__Ronggon_Asif.csv` — complete combined analysis dataset, including original variables and derived labels/log rate.
3. `M5Assign 5.1 Report__Ronggon_Asif.pdf` — dataset description, results and interpretation, figures/tables, and conclusion.

The script and CSV use the required `Data__Last Name_First Name` base; the PDF uses `Report__Last Name_First Name`. The project retains the supplied `Ronggon_Asif` name. The source CSVs and dataset attribution are also preserved in the repository. See `data/README.md` for the data license and original-file checksums.

## Final checks before uploading

- [ ] Confirm the student's name and course details are correct in the report and filenames.
- [ ] Confirm the latest R/LaTeX build succeeded and that the PDF contains the latest results.
- [ ] Open the PDF; inspect plot labels, tables, page breaks, figure references, and the confusion matrix at normal reading size.
- [ ] Confirm both LDA **and** multivariate regression appear in the report and R script.
- [ ] Confirm the `.csv` opens with 2,110 observations and the required quantitative and categorical variables.
- [ ] Rerun the submitted `.R` with the submitted `.csv` in the same working directory; confirm it regenerates the analysis. Preserve `results/sessionInfo.txt` as the environment record.
- [ ] Upload all three required files to Moodle and verify its submission status. GitHub publication does not perform this step.
- [ ] Submit by the course deadline: **Friday, September 25, 2026, 11:59 p.m.** The supplied instructions call the time zone **CST**; check the course's displayed time-zone setting.
