# ENCODE primary human VTRNA1 expression atlas

A census of **VTRNA1-1, VTRNA1-2 and VTRNA1-3** in every released ENCODE human primary-cell and
primary-tissue RNA-seq dataset that has a usable gene-level quantification. It is built only from
ENCODE: no GTEx, cell lines, organoids, tumours or external data.

**Experimentally manipulated primary cells and tissues are included** and labelled. These are
samples treated or stimulated in vitro, activated or polarized, genetically modified, or flagged
`perturbed` by ENCODE. Every row carries `manipulation_status` (`baseline` /
`possibly_manipulated` / `manipulated`) and `manipulation_details`. To get the untreated set,
filter on `manipulation_status == "baseline"`, or run with `--exclude-manipulated`, which drops
them entirely and reproduces the earlier 725-experiment census.

* **ENCODE queried:** 2026-11-25 (UTC), via the REST API at `https://www.encodeproject.org`.
  Each raw response is cached with its retrieval timestamp in `data/encode_metadata/api_cache/`,
  and every request is logged in `data/encode_metadata/query_log.tsv`.
* **Reproduce/update:** `python3 download_encode_vtrna1.py` uses the cache. Add `--refresh` to
  re-query ENCODE.

---

## Headline result: library type decides whether vault RNAs are visible at all

| Track | Assay / pipeline | Bio. samples | VTRNA1-1 detected | VTRNA1-2 detected | VTRNA1-3 detected |
|---|---|---|---|---|---|
| **A** total RNA-seq, ENCODE uniform pipeline (STAR + **RSEM**, GENCODE V29) | Caltech, CSHL, UW, Stanford | 559 (505 baseline, 3 possibly, 51 manipulated) | TPM>0 in **0/559** | 0/559 | 0/559 |
| **B** total RNA-seq, UW/Altius lab pipeline (**STAR unique-read counts**, GENCODE V43) | UW immune primary cells | 375 (285 baseline, 14 possibly, 76 manipulated) | count>0 in 362/375 (median CPM 0.43) | 57/375 (median 0) | 368/375 (median CPM 1.00) |
| **C** small RNA-seq (<200 nt), ENCODE small-RNA pipeline (STAR counts, V24) — *reported separately* | CSHL tissues & primary cells | 129 | 129/129 (median CPM **35.8**) | 114/129 (median 1.9) | 108/129 (median 0.43) |

**The zeros in Track A are an artefact of the method, not evidence that the genes are silent.**
In every ENCODE RSEM `genes.results` file, VTRNA1-1 is 99 nt long and VTRNA1-2/1-3 are 89 nt long.
All three are shorter than the library's mean fragment length, so RSEM sets
`effective_length = 0`. It then forces `expected_count`, `TPM` and `FPKM` to 0 whatever reads
are present. This happens in 557/559 files for VTRNA1-1 and 559/559 for VTRNA1-2/1-3. The other
two files (placenta and thymus, with 120–200 nt inserts) have effective lengths of about 1 and
still give 0. These rows are kept with their reported values. They are labelled
`measurement_status = present_effective_length_0`, and `quantifiable = False` with an explicit
`na_reason`. **Do not interpret Track A TPMs as expression.**

Track C, small RNA-seq of the <200 nt fraction with TAP treatment so 5′-triphosphate Pol III
transcripts can be captured, is the only ENCODE dataset class that is clearly appropriate for
intact ~100-nt vault RNAs. It shows VTRNA1-1 ≫ VTRNA1-2 > VTRNA1-3 across tissues. Track B gives
low but non-zero unique-read counts from total-RNA libraries.

**The tracks measure different things in different units and must not be pooled.**

---

## Output files

| File | Content |
|---|---|
| `encode_primary_human_vtrna1_expression.csv` | **Primary output.** One row per biological replicate (selected gene-quantification file) of total RNA-seq (Tracks A + B; see `quant_schema`). 934 rows, 134 columns, of which 790 are baseline, 17 possibly manipulated and 127 manipulated. |
| `encode_primary_human_vtrna1_expression_long.csv` | Long format of the above: sample × gene (2,802 rows). |
| `encode_primary_human_vtrna1_smallRNA_expression.csv` / `_long.csv` | Track C, kept separate as requested (129 rows / 387 rows). |
| `encode_vtrna1_source_manifest.csv` | **Every experiment considered** (4,052), plus every gene-quant file seen for included experiments. `row_type` is `selected_file`, `alternate_file_not_selected`, `experiment_excluded` or `eligible_no_gene_quantification`, and each non-selected row has an `exclusion_reason`. |
| `vtrna1_gene_identifier_audit.csv` | Per track × assembly × annotation × quantification schema: GTF gene ID, ID as reported in files, stable ID, files present/absent/unparsed, and files with effective length 0. Includes annotations seen only in non-selected files (hg19 V19, GRCh38 V24). |
| `summary_vtrna1_by_gene_*.csv` | Per-gene n measured, n>0, median/mean/max, biosample and file with the maximum. Kept separate for RSEM TPM, UW STAR CPM and small-RNA CPM, and given for the subsets `all`, `baseline`, `possibly_manipulated` and `manipulated`. |
| `summary_vtrna1_by_biosample_*.csv` | Per biosample term × `manipulation_status`: N biological samples, N experiments, N donors, median value and n>0 for each gene, plus the manipulations involved. |
| `qc_report.txt` | Full QC log, also printed by the script. |
| `data/extracted/vtrna1_rows_extracted_from_source_files.csv` | Every VTRNA1 row exactly as parsed from each source file (file accession → raw values). |
| `data/extracted/vtrna1_gtf_records.csv` | VTRNA1 gene records from each ENCODE GENCODE GTF used. |
| `data/encode_metadata/` | Cached API responses (gzipped JSON with URL and timestamp), query log, run info. |
| `data/quantifications/<track>/<assembly>/ENCFF*.tsv.gz` | Raw downloaded quantification files, md5-verified against ENCODE and stored gzipped (not committed to git; about 1.8 GB; rebuilt by the script). |
| `data/references/ENCFF*.gtf.gz` | ENCODE GENCODE GTFs (not committed; rebuilt by the script). |
| `data/logs/` | `run.log`, `failed_downloads.tsv` (empty in this run). |

Key columns in the wide table:

* `VTRNA1_x_gene_id` is the identifier exactly as written in the source TSV.
* `VTRNA1_x_TPM`, `_expected_count`, `_FPKM`, `_length`, `_effective_length` hold RSEM values.
  They are NA for STAR-count files, which have no such values.
* `VTRNA1_x_count_sense`, `_CPM` and `_count_unstranded/_strand1/_strand2` hold STAR values.
  They are NA for RSEM files.
* `VTRNA1_x_measurement_status`, `_quantifiable` and `_na_reason` record measurement validity.
* `quant_schema`, `library_rna_class`, `vault_rna_suitability`, `library_size_range`,
  `library_depleted_in`, `library_treatments`, `library_documents`, `strandedness`, `read_length`,
  `run_type`, `sequencing_platform`, `analysis_title`, `pipeline_version` and `read_depth_qc` give
  technical context.
* `manipulation_status` / `manipulation_details` record experimental manipulation (see Inclusion criteria).
* `experiment_audit*` holds ENCODE audit levels and categories. `experiment_flags` holds disease
  status, possible in-vitro polarization and library-prep enzymes.
* `n_rows_same_donor_and_biosample_term` counts rows sharing the same donor and biosample term
  (see Limitations).

---

## ENCODE API queries

All queries use `limit=all`, and `frame=object` unless `field=` selection is used:

```
/search/?type=Experiment&status=released&replicates.library.biosample.donor.organism.scientific_name=Homo sapiens
         &assay_title=total RNA-seq&assay_title=small RNA-seq&assay_title=polyA plus RNA-seq&assay_title=polyA minus RNA-seq
         &assay_title=microRNA-seq&assay_title=shRNA RNA-seq&assay_title=siRNA RNA-seq&assay_title=CRISPR RNA-seq
         &assay_title=CRISPRi RNA-seq&assay_title=long read RNA-seq
/search/?type=Replicate&experiment.assay_title=<same 10 titles>
/search/?type=Library&nucleic_acid_term_name=RNA            (+ explicit batch fetch of replicate-referenced libraries missing from this query)
/search/?type=Biosample&organism.scientific_name=Homo sapiens
/search/?type=BiosampleType   /search/?type=HumanDonor   /search/?type=Treatment   /search/?type=Platform   /search/?type=Lab
/search/?type=File&output_type=gene quantifications&assay_title=<same 10 titles>
/search/?type=File&file_format=fastq&assay_title=total RNA-seq&assay_title=small RNA-seq&field=accession&field=dataset&field=biological_replicates
         &field=technical_replicates&field=read_length&field=run_type&field=paired_end&field=platform&field=status ...
/search/?type=File&file_format=gtf&output_type=genome reference
/search/?type=Analysis&accession=ENCAN...   (batches of 80; analyses of eligible experiments)
/search/?type=Document&uuid=...             (batches of 60; library protocol documents)
/files/<ENCFF>/@@download/<ENCFF>.tsv       (quantifications; redirects to encode-public.s3.amazonaws.com)
/files/<ENCFF>/@@download/<ENCFF>.gtf.gz    (GENCODE GTFs ENCFF159KBI, ENCFF824ZKD, ENCFF271PVT)
```

---

## Inclusion criteria

An experiment is included when **all** of the following hold:

1. Organism *Homo sapiens*, `status = released`.
2. `assay_title = "total RNA-seq"` (ENCODE `assay_term_name = RNA-seq`; Tracks A/B). Small RNA-seq
   experiments go into the separate Track C.
3. `biosample_ontology.classification` is `primary cell` or `tissue`. This is checked for the
   experiment **and for every replicate's biosample**.
4. Experimental manipulation does **not** exclude a sample by default. It is recorded instead:
   * **`manipulated`**: any of
     * biosample `treatments` (anti-CD3/CD28 beads, IL-2, IL-4, IL-12/15/18, CpG, calcium,
       polarizing cytokines and antibodies)
     * `genetic_modifications` or `applied_modifications`
     * a biosample or experiment `perturbed` flag
     * a biosample term that by definition denotes an in vitro activated, stimulated or polarized
       state ("activated …", "stimulated …", "inflammatory macrophage", "suppressor macrophage").
       ENCODE records no `Treatment` for many UW samples of this kind, so the term itself is the
       evidence.
   * **`possibly_manipulated`**: T-helper, Treg and Tfh subsets with no recorded treatment. They
     may have been sorted ex vivo or polarized in vitro.
   * **`baseline`**: everything else.

   `--exclude-manipulated` turns `manipulated` back into an exclusion reason.
5. Library metadata is consistent with total RNA (or small RNA for Track C). Libraries depleted of
   polyadenylated RNA are rejected from the total track.
6. At least one **released** gene-quantification TSV covering a single biological replicate.

Library metadata (`depleted_in_term_name`, `size_range`, `library_size_selection_method`,
`extraction_method`, `fragmentation_methods`, `strand_specificity`, protocol documents) is
recorded for every row. It is used to classify library type, not just the assay title.

## Exclusion criteria (each recorded in the manifest)

* **Cell lines**, including the 599 Montgomery-lab AFGR lymphoblastoid lines.
  `BiosampleType` objects missing from the bulk search have their classification taken from the
  ontology `@id`, which encodes it.
* **In-vitro differentiated cells and organoids.**
* **Tumour/cancer biosamples.** Disease or term regex; none were found among primary
  cells/tissues.
* **Immortalized, transformed, xenograft or iPSC descriptions** in biosample text, or a `host`
  organism.
* **In-vitro treated, stimulated or polarized biosamples**, **only** when run with
  `--exclude-manipulated`. By default they are included (see above).
* **Other assays**, each with its own reason:
  * poly(A)+ libraries: vault RNAs are not polyadenylated.
  * poly(A)− fraction.
  * microRNA-seq: <30 nt libraries with miRNA-only quantification.
  * shRNA, siRNA, CRISPR and CRISPRi perturbation assays.
  * long-read RNA-seq.
* **Experiments where a replicate's library or biosample cannot be resolved.** This guards
  against silently skipping checks.

**Not excluded:**

* Experimentally manipulated primary cells/tissues: 119 total RNA-seq experiments and 127
  samples are labelled `manipulated`. No small RNA-seq samples are manipulated. The
  knockdown/CRISPR assays remain excluded, because every one of them uses cell lines.

* Disease-annotated non-cancer tissue: dorsolateral prefrontal cortex from donors with
  Alzheimer's disease or (mild) cognitive impairment, and one coronary-artery-disease donor.
  These are listed in `disease_term_name`, `health_status` and `experiment_flags`.
* Library-prep enzymatic "treatments" (Tobacco Acid Pyrophosphatase, Antarctic Phosphatase).
  These are recorded in `library_treatments`. DSN normalization would be flagged, but none occurs
  in the included data.

## Classification of primary cells and tissues

ENCODE's own `BiosampleType.classification` is authoritative. `organ` comes from `organ_slims`,
`system` from `system_slims` and `cell_type` from `cell_slims` (primary cells only).
`anatomical_site` is the tissue term name. Donor sex, age, age units and life stage come from the
biosample's own fields, falling back to the donor's. `donor_ethnicity` is reported only where the
ENCODE `HumanDonor.ethnicity` field is filled in. Nothing is inferred.

## Gene-ID resolution

The quantification TSVs contain **no gene symbols**. Symbols were therefore resolved from the
GENCODE GTF that ENCODE distributes for each annotation version:

| Assembly / annotation | ENCODE GTF | VTRNA1-1 | VTRNA1-2 | VTRNA1-3 |
|---|---|---|---|---|
| GRCh38 V29 (ENCODE4 RSEM) | ENCFF159KBI | ENSG00000199990.1 | ENSG00000202111.1 | ENSG00000202515.1 |
| GRCh38 V24 (small RNA; ENCODE3 alternates) | ENCFF824ZKD | ENSG00000199990.1 | ENSG00000202111.1 | ENSG00000202515.1 |
| hg19 V19 (alternates only) | ENCFF271PVT | ENSG00000199990.1 | ENSG00000202111.1 | ENSG00000202515.1 |
| GRCh38 **V43** (UW lab pipeline; version from file `submitter_comment`) | *none on ENCODE* | matched by stable ID | matched by stable ID | matched by stable ID |

The GTFs were scanned for `gene_name` starting `VTRNA1-` and for the seed ID ENSG00000199990.
Each quantification file was then matched by **version-stripped Ensembl ID**, using the IDs from
its own annotation plus the union across annotations, and by a symbol column where one exists
(none do). The original identifier is preserved in `VTRNA1_x_gene_id`. The identifiers are
**identical across V19, V24 and V29**, and every one of the 1,063 parsed files contains all three
rows. Apparent absence is therefore never caused by an annotation change.

For V43, ENCODE hosts no GTF. Gene symbols for those 375 files rely on the stable IDs being
unchanged, which is documented in `VTRNA1_x_gene_symbol_source` and in the audit.

## Biological vs technical replicates; choosing among multiple files

* **One row = one biological replicate = one gene-quantification file.**
* The analysis used is the experiment's ENCODE **`default_analysis`**. If that has no released
  gene quantification, the order of preference is GRCh38, then newest GENCODE, then the uniform
  pipeline, then the newest file. In this run, all selected files come from the default analysis
  where one exists.
* Within the chosen analysis, files covering more than one biological replicate (pooled) are not
  used. If one biological replicate has several files, ENCODE's **`preferred_default`** file is
  kept, then the newest. This happens in 16 rows, where a biological replicate has several
  technical-replicate libraries each quantified separately (for example ENCSR584CVV, with 10
  technical replicates). The unused files are listed in `other_technical_replicate_files_not_selected`
  and in the manifest.
* hg19 V19 / GRCh38 V24 outputs, reruns, archived and revoked files are all listed as
  `alternate_file_not_selected` with a reason (606 total-RNA and 129 small-RNA alternate files).
* Library, fastq, read-length and platform metadata are restricted to the technical replicate(s)
  in the selected file.
* No values are averaged.

## Zero vs NA

| Situation | Value | `measurement_status` |
|---|---|---|
| Gene row present, value 0, gene quantifiable | `0` | `present` |
| Gene row present, RSEM `effective_length = 0` | reported `0` kept | `present_effective_length_0`, `quantifiable = False`: **structural zero, not a measurement** |
| Metric not produced by that pipeline (e.g. TPM in STAR-count files) | `NA` | `present` (the counts are the measurement) |
| Gene absent from file although in the annotation | `NA` | `NA_absent_from_file` (never occurred) |
| Gene absent from annotation and file | `NA` | `NA_absent_from_annotation_and_file` (never occurred) |
| Download or parse failure | `NA` | `NA_file_unavailable` / `NA_parse_failed` (never occurred) |

**CPM** is derived by this pipeline: sense-strand count divided by the total sense-strand counts
over all genes in that file, × 10⁶. The sense column is chosen empirically per file: the strand
column whose genome-wide total is more than 3× the other, otherwise the unstranded column. It came
out reverse for the UW files and forward for CSHL small RNA, matching library strandedness.

---

## QC summary (from `qc_report.txt`)

| | Total RNA-seq (A+B) | Small RNA-seq (C) |
|---|---|---|
| ENCODE experiments queried (all 10 assay titles) | 4,052 | — |
| Experiments with this assay | 1,650 | 176 |
| Passing inclusion criteria | **844** (708 baseline, 17 possibly manipulated, 119 manipulated) | **91** (all baseline) |
| Eligible but no usable quantification | 0 | 0 |
| Unique biosample terms | 172 | 57 |
| Unique donors | 298 | 53 |
| Biological samples / files selected | 934 (790 / 17 / 127) | 129 |
| Downloaded (md5 OK) / parsed | 934 / 934 | 129 / 129 |
| Files containing VTRNA1-1 / -2 / -3 rows | 934 / 934 / 934 | 129 / 129 / 129 |
| Value > 0, VTRNA1-1 | TPM: 0/559; STAR count: 362/375 | 129/129 |
| Value > 0, VTRNA1-2 | TPM: 0/559; STAR count: 57/375 | 114/129 |
| Value > 0, VTRNA1-3 | TPM: 0/559; STAR count: 368/375 | 108/129 |
| Annotations | GRCh38 V29 (RSEM, 559); GRCh38 V43 (STAR, 375) | GRCh38 V24 (129) |
| Files with any gene missing | 0 | 0 |
| Duplicate file / (experiment, bio-rep) / library / biosample rows | 0 / 0 / 0 / 0 | 0 / 0 / 0 / 0 |
| Failed downloads | 0 | 0 |

Parsing sanity checks:

* Sum of TPM per RSEM file was between 999,994 and 1,000,001.
* Every RSEM file has 59,526 genes.
* No negative values.
* Spot checks against the raw TSVs matched.

### Per-gene summary

*Track A (RSEM TPM, n = 559; same for every manipulation subset):* every gene has median, mean and maximum TPM of 0, because the
values are structural zeros (see above).

*Track B (UW total RNA-seq, STAR unique-read CPM, immune primary cells):*

| Subset | Gene | n | n>0 | median | mean | max | max biosample |
|---|---|---|---|---|---|---|---|
| all | VTRNA1-1 | 375 | 362 | 0.43 | 0.68 | 5.09 | CD14-positive monocyte |
| all | VTRNA1-2 | 375 | 57 | 0.00 | 0.013 | 0.49 | suppressor macrophage |
| all | VTRNA1-3 | 375 | 368 | 1.00 | 1.55 | 11.74 | CD14-positive monocyte |
| baseline | VTRNA1-1 | 285 | 275 | 0.39 | 0.65 | 5.09 | CD14-positive monocyte |
| baseline | VTRNA1-2 | 285 | 29 | 0.00 | 0.007 | 0.27 | CD14-positive monocyte |
| baseline | VTRNA1-3 | 285 | 280 | 1.18 | 1.65 | 11.74 | CD14-positive monocyte |
| manipulated | VTRNA1-1 | 76 | 75 | 0.59 | — | 3.89 | suppressor macrophage |
| manipulated | VTRNA1-2 | 76 | 27 | 0.00 | — | 0.49 | suppressor macrophage |
| manipulated | VTRNA1-3 | 76 | 74 | 0.66 | — | 5.75 | T-helper 1 cell |

In vitro polarized macrophages (inflammatory n = 18, suppressor n = 17) have the highest
VTRNA1-1 (median CPM 0.98 and 1.17) and VTRNA1-2 of any manipulated group, but lower VTRNA1-3
(0.42 and 0.37) than T cells.

*Track C (small RNA-seq, CPM, n = 129):*

| Gene | n | n>0 | median | mean | max | max biosample |
|---|---|---|---|---|---|---|
| VTRNA1-1 | 129 | 129 | 35.8 | 50.7 | 231.4 | right lobe of liver |
| VTRNA1-2 | 129 | 114 | 1.88 | 16.4 | 144.6 | placental pericyte |
| VTRNA1-3 | 129 | 108 | 0.43 | 1.67 | 54.6 | placental pericyte |

Tissue and cell-type medians for every biosample term are in `summary_vtrna1_by_biosample_*.csv`
and at the end of `qc_report.txt`.

---

## Known limitations

1. **RSEM structurally cannot quantify these loci.** The ENCODE uniform long-RNA pipeline
   (Track A, 508 samples, most tissues) gives no usable VTRNA1 information. A real measurement
   would require re-quantifying the BAMs, for example with featureCounts on unique alignments,
   which was out of scope ("use ENCODE quantification files").
2. **Library size selection.**
   * 248 total-RNA libraries (CSHL) are explicitly size-selected >200 nt, which largely excludes
     mature vault RNAs.
   * 549 (Caltech, made with the NEBNext Ultra II Directional kit, and UW) record no depletion or size
     selection. Their `vault_rna_suitability` is "uncertain".
   * Even without explicit selection, bead cleanups and fragment-size distributions deplete
     ~100-nt RNAs. Total-RNA counts are therefore lower bounds and not comparable across labs.
3. **Paralog mappability.** VTRNA1-1/-2/-3 are highly similar paralogs in one cluster at
   5q31.3. Both STAR-count pipelines count only uniquely mapped reads, so per-paralog counts
   reflect uniquely mappable positions. Multimapping reads, reported only in aggregate as
   `N_multimapping`, are discarded. The paralog ratios differ between Track B
   (VTRNA1-3 > VTRNA1-1) and Track C (VTRNA1-1 ≫ VTRNA1-2 > VTRNA1-3), probably for this reason
   and because of library chemistry.
4. **Track B is immune cells only**, all from one lab, and is concentrated in a few donors.
   * 76 donor × cell-type groups have several rows (441 rows in total). ENCDO116ULP, ENCDO833KUM and ENCDO923WRX
     contribute 38–47 CD4⁺ T-cell samples each.
   * These are distinct ENCODE biosamples, but ENCODE records no condition or time point that
     separates them.
   * They are not independent donors. Use `donor_accession` and
     `n_rows_same_donor_and_biosample_term` for donor-level analyses.
5. **UW metadata gaps.**
   * UW libraries are recorded with `nucleic_acid_term_name = DNA`, flagged in
     `experiment_flags`.
   * They have no strand, size or depletion fields.
   * The annotation version (V43) comes from free-text `submitter_comment`.
   * T-helper subsets may be polarized in vitro (flagged).
6. **Small RNA-seq (Track C)** is ENCODE2/3 CSHL data processed only against GRCh38 V24 and
   hg19. Its counts are not length-normalized, and its CPM is relative to the small-RNA
   composition of each library: miRNA, snoRNA and tRNA dominate and differ between tissues.
7. **Manipulated samples are not matched to their own baselines.**
   * Same-donor baseline rows exist for all 47 Caltech and 65 of 76 UW manipulated samples
     (link them by `donor_accession`). But for many UW activated or polarized samples, ENCODE
     gives no treatment, dose or time.
   * Treat `manipulated` rows as additional cell states, not as controlled perturbation
     contrasts.
   * Nearly all manipulated samples are immune cells; the other 4 are calcium-treated
     keratinocytes.
8. **The Alzheimer's/MCI cortex samples** are included as primary tissue. Filter them on
   `disease_term_name` if only non-disease tissue is wanted.
9. **The classification depends on ENCODE metadata as of the query date.** Rerun with
   `--refresh` to pick up new releases. `query_log.tsv` and `run_info.json` record each run.
10. **VTRNA2-1 (nc886) is deliberately not reported.**
