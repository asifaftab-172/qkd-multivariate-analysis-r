# Dataset provenance

Source: Horgan et al. (2026), *QKD Network Characterisation and Coexistence Datasets with Standardised Performance Metrics and Experimental Methods*, v2. [Zenodo record and citation](https://zenodo.org/records/21792844), DOI [10.5281/zenodo.21792844](https://doi.org/10.5281/zenodo.21792844). Published August 4, 2026.

The record lists Jerry Horgan, Niall Quigley, Bennouri Hajar (encoded as "Hajar, Bennouri"), Dmitrii Briantcev, Aleksandra Kaszubowska-Anandarajah, and Daniel Kilper as creators. Consult the source record for authoritative citation metadata.

These five original CSVs contain the 0-km coexistence measurements used by the R script. Their bytes are unchanged. Every downloaded file was checked against the MD5 checksum published in [Zenodo's record metadata](https://zenodo.org/api/records/21792844).

| Condition | Source file | Observations | Source MD5 |
| --- | --- | ---: | --- |
| Baseline | [baseline_no_signal_no_added_fibre_full_sweep.csv](https://zenodo.org/records/21792844/files/baseline_no_signal_no_added_fibre_full_sweep.csv?download=1) | 471 | `df5e9700d992a6807a181d560967198e` |
| 3 dBm | [roadm_3dbm_injected_signal_no_added_fibre_full_sweep.csv](https://zenodo.org/records/21792844/files/roadm_3dbm_injected_signal_no_added_fibre_full_sweep.csv?download=1) | 400 | `e400af6ae4fe693e42a160f9429a8895` |
| 7 dBm | [roadm_7dbm_injected_signal_no_added_fibre_full_sweep.csv](https://zenodo.org/records/21792844/files/roadm_7dbm_injected_signal_no_added_fibre_full_sweep.csv?download=1) | 354 | `6d12509c4306eda78248c21b73951810` |
| 9 dBm | [roadm_9dbm_injected_signal_no_added_fibre_full_sweep.csv](https://zenodo.org/records/21792844/files/roadm_9dbm_injected_signal_no_added_fibre_full_sweep.csv?download=1) | 319 | `847d94e867b7e6ec96105a52c3025d88` |
| 12 dBm | [roadm_12dbm_injected_signal_no_added_fibre_full_sweep.csv](https://zenodo.org/records/21792844/files/roadm_12dbm_injected_signal_no_added_fibre_full_sweep.csv?download=1) | 566 | `9a7b916047e0d04caeae34567ef62a33` |

Total: **2,110 observations**. Original columns are `Noise_dBm`, `Sweep_km`, `VOA_dB`, `QBER`, and `SecureKeyRate_bps`. QBER is a fraction, and secure key rate is measured in bits per second. The analysis adds condition labels and `log10(SecureKeyRate_bps + 1)` in the derived combined file under `data/processed/`.

## License

The source data are licensed by their creators under [Creative Commons Attribution 4.0 International (CC BY 4.0)](https://creativecommons.org/licenses/by/4.0/). Retain attribution, the source DOI, and this license link when redistributing. The repository's code license does not replace the dataset license.
