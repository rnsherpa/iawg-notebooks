# HepG2 MPRA versus liver eQTL

Zhiyuan (Andy) Chu (`zychu912`). Two notebooks compare Mohlke and Tewhey HepG2 MPRA effects with GTEx liver eQTLs, including FINEMAP/SuSiE and top/all eGene comparisons.

## Files and running

- `notebooks/`: two executed notebooks and one shared helper, `eqtl_mpra.py`.
- `plots/`: six PNG figures at 300 DPI.
- `requirements.txt`: tested package versions. The recorded run used Python 3.9.6.

Install the dependencies in your Python environment with `python -m pip install -r requirements.txt`. Set `IAWG_DATA_DIR` to your input directory, open a notebook with that environment in Jupyter, and run all cells. If needed, register the environment with `python -m ipykernel install --user --name iawg --display-name IAWG` and select that kernel. The notebooks contain the summary tables, input checksums and run information; separate result files are unnecessary. Plot files are overwritten when rerun; inputs are never modified.

## Required inputs

| File | Source |
|---|---|
| `Mohlke_HepG2_MPRA_LipidTraits.tsv` | [Synapse syn68541428.2](https://www.synapse.org/Synapse:syn68541428.2) |
| `GTEx_Liver_eQTLs.signif_pairs.tsv.gz` | [Synapse syn69058503.1](https://www.synapse.org/Synapse:syn69058503.1) |
| `gtex_finemap_eqtls_pip01_with_spdi.csv` | Existing processed table; original download and SPDI-conversion provenance remain to be documented. |
| `hepg2_combined_final_processed.tsv.gz` | Existing processed Tewhey table; processing provenance and the supplied emVar definition remain to be documented. |

Raw data are not included. Saved outputs can be viewed without downloading inputs. Full re-execution requires these files; this is not yet an end-to-end public download pipeline.

## Methods and limitations

- GTEx: nominal p below the supplied gene threshold. Fine-mapping: Liver, PIP>0.1, using marginal beta.
- Top eGene: minimum nominal p (GTEx) or maximum PIP (fine-mapping) per variant; ties resolve by gene ID, effect, then source row ID.
- Highlighting: Mohlke q<0.05 or supplied Tewhey emVar. The same mask defines plot colors and subgroup statistics.
- Distinct MPRA records are retained. Rows are record × association matches; unique variants, variant–gene pairs and MPRA records are counted separately. Correlations are descriptive, without independent-sample p-values.
- Concordance requires nonzero effects of the same sign; zero-effect rows remain in the denominator. Invalid labels fail explicitly; nonfinite effect/score exclusions are counted.
- Effects are reported as supplied, **not independently direction-harmonized**. Some fine-mapping allele fields disagree with hg38 SPDI. The notebooks include a metadata-consistent sensitivity subset; literal allele agreement does not establish the effect-column convention, and unresolved indels may be excluded from that subset. No signs are automatically flipped.

These comparisons do not establish causality or evaluate a variant prediction model. Source direction definitions and incomplete provenance require follow-up. Validation: both notebooks execute from clean kernels, and seven local regression tests cover data handling. The validation harness is retained locally rather than included in this contribution.
