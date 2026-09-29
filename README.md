# CDC14B circRNA-miRNA axis in breast cancer

A reproducible BIOL 260 bioinformatics project that reconstructs a **hypothesis-generating regulatory axis** from the assigned cell-cycle gene **CDC14B** to a CDC14B-derived circular RNA, a breast-cancer-relevant miRNA, and predicted downstream mRNA targets.

## Why this project exists

The original take-home exercise asks for six steps:

1. start from an assigned gene;
2. identify two gene-derived circular RNAs;
3. predict three miRNAs that could interact with one of those circRNAs;
4. choose one miRNA and examine its breast-cancer relevance and survival association;
5. predict three downstream mRNA targets with TargetScan;
6. summarize the resulting axis.

Rather than treating those steps as independent database lookups, this repository connects them into one explicit biological question:

> **Could a CDC14B-derived circRNA plausibly alter the availability of a breast-cancer-relevant miRNA and thereby influence repression of downstream mRNA targets?**

The answer is presented as a **computational hypothesis**, not as proof of a physical miRNA-sponge mechanism.

## Final proposed axis

```text
CDC14B
  |
  v
hsa_circ_0087640
  |
  | putative miRNA sequestration
  v
miR-136-5p
  |
  | predicted miRNA-mediated repression
  v
ORC6 / XRCC2 / TFAP2C
```

Compact notation:

```text
CDC14B -> hsa_circ_0087640 -| miR-136-5p -| {ORC6, XRCC2, TFAP2C}
```

The two inhibitory edges are **predicted regulatory relationships**, not experimentally demonstrated interactions in this project.

## Biological idea in plain language

A **circular RNA (circRNA)** is produced by back-splicing, which joins a downstream splice donor to an upstream splice acceptor and creates a covalently closed RNA molecule. Some circRNAs contain sequences complementary to the **seed region** of a microRNA (miRNA).

A **miRNA** is a short regulatory RNA. Its nucleotides 2-8 form the seed, a major determinant of target recognition. When a miRNA binds a compatible site on an mRNA, it can reduce the expression of that mRNA.

If a circRNA contains accessible miRNA-binding sites and is present at sufficient abundance in the correct cellular compartment, it may compete for that miRNA. This is often described as **miRNA sponging**. Importantly, sequence compatibility alone does not prove sponging; abundance, localization, accessibility, and stoichiometry matter.

## Why the workflow is more than a database lookup

The linked CircInteractome page in the assignment was unavailable during the analysis. Instead of replacing it with an unrelated predictor, the repository reconstructs the intended sequence-based logic and makes every step inspectable.

The analysis follows a **sequence-first evidence-integration workflow**:

```text
CDC14B
  -> circRNA annotation
  -> canonical miRNA seed-site scan
  -> seed-family collapse
  -> breast-cancer evidence
  -> Kaplan-Meier analysis
  -> TargetScan downstream targets
  -> proposed regulatory axis
```

The extra work is mainly a quality-control layer: it reduces cherry-picking, separates sequence strength from disease relevance, and makes the limitations visible.

## Main result

Three miRNAs were retained as the principal candidates:

| Candidate | Sequence evidence | Why it remained relevant |
|---|---:|---|
| `hsa-miR-217-5p` | strongest of the three (`S_seed = 7`) | independent breast-cancer tumor-suppressive evidence |
| `hsa-miR-106a-3p` | `S_seed = 5` | mammary epithelial plasticity / tumoroid-initiation evidence |
| `hsa-miR-136-5p` | `S_seed = 5` | breast-cancer relevance, circRNA-sponging precedent, and METABRIC availability |

`miR-136-5p` was used for the formal homework axis because it could be followed through **all required layers**, including the assignment-specified METABRIC survival analysis. This does **not** mean it had the strongest sequence score; `miR-217-5p` was the strongest sequence-derived candidate.

### Kaplan-Meier result

KM Plotter settings:

- miRNA identifier: `hsa-miR-136`
- cohort: METABRIC
- endpoint: overall survival
- split: median
- follow-up: unrestricted
- clinical restrictions: none

Result:

```text
N = 1262
low expression  = 637
high expression = 625
HR = 0.84
95% CI = 0.69-1.02
log-rank p = 0.081
```

The point estimate is directionally favorable for higher miR-136 expression, but the result is **not statistically significant** because the confidence interval includes 1 and `p > 0.05`.

## Repository layout

```text
.
├── BIOL260_CDC14B_report.pdf
├── BIOL260_CDC14B_report.tex
├── references.bib
├── README.md
├── .gitignore
├── data/
│   ├── cdc14b_circ_0087640.fasta
│   ├── cdc14b_circrna_candidates.tsv
│   ├── sequence_shortlist.csv
│   ├── 04_hsa-miR-136-5p_TargetScan.csv
│   ├── targetscan_gene_level.csv
│   ├── candidate_target_summary.csv
│   ├── cancer_theme_summary.csv
│   └── km_mir136_metabric.txt
├── figures/
│   └── km_mir136_metabric.pdf
└── scripts/
    ├── 01_scan_cdc14b_circrna_mirna.py
    ├── 02_collapse_seed_families.py
    ├── 03_query_targetscan_multimir.R
    └── 04_target_network_go.R
```

## Reproducing the analysis

### Requirements

The Python scripts use only the Python standard library.

The R scripts use Bioconductor packages. The TargetScan query uses `multiMiR`; the network/GO analysis uses annotation packages such as `AnnotationDbi`, `org.Hs.eg.db`, and `GO.db`.

### 1. Scan the circRNA against mature human miRNAs

```bash
python scripts/01_scan_cdc14b_circrna_mirna.py \
  --circ-fasta data/cdc14b_circ_0087640.fasta \
  --out stage1/cdc14b_circ_0087640_mirna_scan.csv \
  --top 0
```

If `--mirna-fasta` is omitted, the script downloads the pinned human mature-miRNA reference used in the project.

The scan searches canonical 8mer, 7mer-m8, and 7mer-A1 sites and also searches across the **back-splice junction**.

### 2. Collapse mature miRNAs into seed families

```bash
python scripts/02_collapse_seed_families.py \
  --summary stage1/cdc14b_circ_0087640_mirna_scan.csv \
  --sites stage1/cdc14b_circ_0087640_mirna_scan_sites.csv \
  --outdir stage2
```

The ranking heuristic is:

```text
S_seed = 3*(8mer) + 2*(7mer-m8) + 1*(7mer-A1)
```

It is a transparent prioritization score, **not** a validated biophysical binding score.

### 3. Query TargetScan through multiMiR

```bash
Rscript scripts/03_query_targetscan_multimir.R stage3_targetscan
```

The script queries the three final candidates, maps `miR-217-5p` to the older TargetScan alias `hsa-miR-217`, and keeps the pre-specified top 20% of TargetScan predictions.

### 4. Gene-level consolidation and exploratory GO analysis

```bash
Rscript scripts/04_target_network_go.R \
  stage3_targetscan/targetscan_finalists_combined.csv \
  stage4
```

`04_target_network_go.R` collapses TargetScan site/transcript records to gene-level predictions and performs the exploratory background-sensitivity analysis used in the report.

The GO analysis is **exploratory**. Apparent enrichment under an all-human-gene background was not robust to a TargetScan-union background, so no pathway-level mechanism was used to justify the final axis.

### 5. Kaplan-Meier analysis

The survival step was performed manually in KM Plotter because the course exercise explicitly requires that framework.

The exported patient-level data are stored in:

```text
data/km_mir136_metabric.txt
```

and the original plot is stored in:

```text
figures/km_mir136_metabric.pdf
```

## How to interpret the final model

The repository supports the following statements:

- `hsa_circ_0087640` is a CDC14B-derived circRNA with a mature sequence suitable for sequence-based screening;
- its sequence contains canonical sites compatible with `miR-136-5p`;
- `miR-136-5p` is independently relevant to breast-cancer biology;
- TargetScan predicts `ORC6`, `XRCC2`, and `TFAP2C` among its potential mRNA targets;
- higher miR-136 expression shows a directionally favorable, but non-significant, overall-survival trend in METABRIC.

It **does not** establish that:

- `hsa_circ_0087640` physically binds miR-136-5p in breast-cancer cells;
- the circRNA functions as an effective sponge in vivo;
- the three selected mRNAs are functionally derepressed by the circRNA;
- miR-136 is a statistically validated prognostic biomarker in METABRIC.

## What would validate the hypothesis experimentally?

A reasonable experimental follow-up would include:

1. junction-specific RT-qPCR to confirm `hsa_circ_0087640`;
2. RNase R resistance and/or sequencing across the back-splice junction;
3. subcellular localization of the circRNA;
4. AGO2-RIP, RNA pull-down, or luciferase assays for circRNA-miRNA binding;
5. perturbation of circRNA/miRNA levels followed by measurement of `ORC6`, `XRCC2`, and `TFAP2C`;
6. rescue experiments to test whether the proposed regulatory relationships are causal.

## Main caveats

- The legacy CircInteractome URL supplied in the assignment returned a page-not-found error, so the sequence-prediction step was reconstructed programmatically.
- The legacy miRCancer interface was not reliably accessible; dysregulation direction was checked against primary breast-cancer literature instead.
- Among the three final candidates, only `hsa-miR-136` was available in the assignment-specified METABRIC miRNA interface at the time of analysis.
- The METABRIC survival association was not statistically significant (`p = 0.081`).
- Sequence and TargetScan predictions are hypothesis-generating, not experimental validation.

## Report

The full explanation, citations, figure, interpretation, and limitations are in:

- [`BIOL260_CDC14B_report.pdf`](BIOL260_CDC14B_report.pdf)
- [`BIOL260_CDC14B_report.tex`](BIOL260_CDC14B_report.tex)

## Citation

If reusing the workflow or code outside the course context, cite the software/databases used in the report, particularly CircInteractome, circBase, TransCirc, TargetScan, multiMiR, and KM Plotter/miRpower. Full references are listed in the report and in `references.bib`.

## Author

**Jad Sawan**  
BIOL 260, American University of Beirut  
September 2026
