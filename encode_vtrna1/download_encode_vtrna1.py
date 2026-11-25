#!/usr/bin/env python3
"""
ENCODE VTRNA1 expression atlas for primary human cells and tissues.

Queries the ENCODE REST API, classifies every released human bulk short-read
RNA-seq-family experiment, selects one gene-quantification TSV per biological
replicate for qualifying primary-cell / tissue experiments, downloads them,
resolves VTRNA1-1 / VTRNA1-2 / VTRNA1-3 gene identifiers from the GENCODE GTF
that ENCODE used for each annotation version, extracts the measurements and
writes the output tables, QC report and summaries.

Two measurement tracks are kept strictly separate:
  * total_RNA  : assay_title "total RNA-seq" (rRNA-depleted, long/total RNA
                 libraries; ENCODE long-RNA pipeline, STAR + RSEM)
  * small_RNA  : assay_title "small RNA-seq" (<200 nt fraction; ENCODE small
                 RNA pipeline, STAR per-gene read counts, no TPM)

Usage:
    python3 download_encode_vtrna1.py                # full run (uses caches)
    python3 download_encode_vtrna1.py --refresh      # re-query ENCODE API
    python3 download_encode_vtrna1.py --workers 8 --outdir .

Re-running is cheap: API responses are cached under data/encode_metadata/,
quantification files are not re-downloaded when a verified copy exists.
"""

from __future__ import annotations

import argparse
import datetime as dt
import gzip
import hashlib
import io
import json
import logging
import os
import re
import sys
import threading
import time
from collections import Counter, defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed
from urllib.parse import urlencode

import numpy as np
import pandas as pd
import requests
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry

BASE = "https://www.encodeproject.org"
HUMAN = "Homo sapiens"

TARGET_SYMBOLS = ["VTRNA1-1", "VTRNA1-2", "VTRNA1-3"]
# Only seed identifier supplied by the user; the others are resolved from GTFs.
SEED_STABLE_IDS = {"VTRNA1-1": "ENSG00000199990"}

PRIMARY_TRACK_ASSAY = "total RNA-seq"
SMALL_TRACK_ASSAY = "small RNA-seq"

# Every assay_title of the bulk/short-read RNA-seq family is "considered" so
# that the manifest documents why each experiment was or was not used.
ASSAY_TITLES_CONSIDERED = [
    "total RNA-seq",
    "small RNA-seq",
    "polyA plus RNA-seq",
    "polyA minus RNA-seq",
    "microRNA-seq",
    "shRNA RNA-seq",
    "siRNA RNA-seq",
    "CRISPR RNA-seq",
    "CRISPRi RNA-seq",
    "long read RNA-seq",
]

ELIGIBLE_CLASSIFICATIONS = {"primary cell", "tissue"}

ASSAY_EXCLUSION_REASON = {
    "polyA plus RNA-seq": "poly(A)+ selected library; vault RNAs are non-polyadenylated Pol III transcripts",
    "polyA minus RNA-seq": "poly(A)-depleted fraction only (not total RNA)",
    "microRNA-seq": "miRNA-size library (<30 nt / 19-25 nt) with miRNA-only quantification; cannot represent intact 88-100 nt vault RNAs",
    "shRNA RNA-seq": "perturbation assay (shRNA knockdown)",
    "siRNA RNA-seq": "perturbation assay (siRNA knockdown)",
    "CRISPR RNA-seq": "perturbation assay (CRISPR genome editing)",
    "CRISPRi RNA-seq": "perturbation assay (CRISPR interference)",
    "long read RNA-seq": "long-read (full-length cDNA) assay, not short-read total RNA-seq gene quantification",
}

CANCER_PATTERN = re.compile(
    r"cancer|carcinoma|tumou?r|neoplas|sarcoma|lymphoma|leuka?emia|glioma|"
    r"glioblastoma|melanoma|blastoma|adenoma|myeloma|metasta|malignan",
    re.I,
)
IMMORTAL_PATTERN = re.compile(
    r"immortali[sz]|transformed|hTERT|SV40|EBV[- ]transformed|"
    r"lymphoblastoid|xenograft|organoid|iPSC|induced pluripotent",
    re.I,
)

# Biosample terms that by definition describe an in vitro activated / stimulated /
# polarized state.  Several such samples carry no Treatment objects in ENCODE, so the
# term itself is used as evidence of experimental manipulation (excluded).
IN_VITRO_STATE_PATTERN = re.compile(r"\bactivated\b|\bstimulated\b|inflammatory macrophage|suppressor macrophage", re.I)
# Terms that may or may not reflect in vitro polarization (flagged, kept).
POSSIBLE_POLARIZED_PATTERN = re.compile(r"T-helper \d+ cell|T-helper cell|regulatory T cell|T follicular helper", re.I)

# ---------------------------------------------------------------------------
# Logging / paths
# ---------------------------------------------------------------------------

log = logging.getLogger("encode_vtrna1")


def utcnow() -> str:
    return dt.datetime.now(dt.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


class Paths:
    def __init__(self, outdir: str):
        self.out = os.path.abspath(outdir)
        self.data = os.path.join(self.out, "data")
        self.meta = os.path.join(self.data, "encode_metadata")
        self.api_cache = os.path.join(self.meta, "api_cache")
        self.refs = os.path.join(self.data, "references")
        self.quant = os.path.join(self.data, "quantifications")
        self.logs = os.path.join(self.data, "logs")
        self.extract = os.path.join(self.data, "extracted")
        for p in [self.api_cache, self.refs, self.quant, self.logs, self.extract]:
            os.makedirs(p, exist_ok=True)


# ---------------------------------------------------------------------------
# ENCODE client with caching + retries
# ---------------------------------------------------------------------------


class EncodeClient:
    def __init__(self, paths: Paths, refresh: bool = False):
        self.paths = paths
        self.refresh = refresh
        self.session = requests.Session()
        retry = Retry(
            total=6,
            backoff_factor=1.5,
            status_forcelist=(429, 500, 502, 503, 504),
            allowed_methods=("GET", "HEAD"),
            raise_on_status=False,
        )
        adapter = HTTPAdapter(max_retries=retry, pool_connections=16, pool_maxsize=16)
        self.session.mount("https://", adapter)
        self.session.headers.update({"Accept": "application/json",
                                     "User-Agent": "encode-vtrna1-atlas/1.0 (python-requests)"})
        self.query_log_path = os.path.join(paths.meta, "query_log.tsv")
        self._lock = threading.Lock()

    def _log_query(self, url: str, cache_file: str, from_cache: bool, n: int | None):
        with self._lock:
            new = not os.path.exists(self.query_log_path)
            with open(self.query_log_path, "a") as fh:
                if new:
                    fh.write("timestamp_utc\turl\tcache_file\tfrom_cache\tn_results\n")
                fh.write(f"{utcnow()}\t{url}\t{os.path.relpath(cache_file, self.paths.out)}\t{from_cache}\t{n}\n")

    def get_json(self, path: str, params: list[tuple[str, str]] | None, cache_name: str) -> dict:
        """GET an ENCODE JSON endpoint, caching the raw response to disk."""
        url = BASE + path + ("?" + urlencode(params) if params else "")
        cache_file = os.path.join(self.paths.api_cache, cache_name + ".json.gz")
        if os.path.exists(cache_file) and not self.refresh:
            with gzip.open(cache_file, "rt") as fh:
                wrapped = json.load(fh)
            self._log_query(url, cache_file, True, len(wrapped["response"].get("@graph", [])) if isinstance(wrapped["response"], dict) else None)
            return wrapped["response"]
        for attempt in range(5):
            try:
                r = self.session.get(url, timeout=600)
                if r.status_code == 404 and "/search/" in path:
                    # ENCODE returns 404 for an empty search result
                    data = r.json()
                    break
                r.raise_for_status()
                data = r.json()
                break
            except (requests.RequestException, ValueError) as exc:
                wait = 5 * (attempt + 1)
                log.warning("GET %s failed (%s); retry in %ss", url, exc, wait)
                time.sleep(wait)
        else:
            raise RuntimeError(f"Failed to GET {url}")
        wrapped = {"url": url, "retrieved_utc": utcnow(), "response": data}
        tmp = cache_file + ".tmp"
        with gzip.open(tmp, "wt") as fh:
            json.dump(wrapped, fh)
        os.replace(tmp, cache_file)
        self._log_query(url, cache_file, False, len(data.get("@graph", [])) if isinstance(data, dict) else None)
        return data

    def search(self, params: list[tuple[str, str]], cache_name: str) -> list[dict]:
        p = list(params) + [("limit", "all"), ("format", "json")]
        data = self.get_json("/search/", p, cache_name)
        return data.get("@graph", [])

    def download(self, href: str, dest_gz: str, md5: str | None) -> tuple[bool, str]:
        """Download a file, verify md5 of the uncompressed bytes, store gzipped.

        Returns (ok, message).  Existing verified files are not re-downloaded.
        """
        sidecar = dest_gz + ".md5"
        if os.path.exists(dest_gz) and os.path.exists(sidecar):
            with open(sidecar) as fh:
                if md5 is None or fh.read().strip() == md5:
                    return True, "cached"
        url = BASE + href
        last_err = ""
        for attempt in range(5):
            try:
                with self.session.get(url, timeout=600, stream=True, allow_redirects=True,
                                      headers={"Accept": "*/*"}) as r:
                    r.raise_for_status()
                    raw = r.content
                # ENCODE md5sum refers to the file exactly as served
                got = hashlib.md5(raw).hexdigest()
                if md5 and got != md5:
                    last_err = f"md5 mismatch (expected {md5}, got {got})"
                    log.warning("%s: %s", url, last_err)
                    time.sleep(3 * (attempt + 1))
                    continue
                os.makedirs(os.path.dirname(dest_gz), exist_ok=True)
                tmp = dest_gz + ".tmp"
                if raw[:2] == b"\x1f\x8b":      # already gzipped: store verbatim
                    with open(tmp, "wb") as fh:
                        fh.write(raw)
                else:                            # plain TSV: store gzipped
                    with gzip.open(tmp, "wb", compresslevel=6) as fh:
                        fh.write(raw)
                os.replace(tmp, dest_gz)
                with open(sidecar, "w") as fh:
                    fh.write(got)
                return True, "downloaded"
            except requests.RequestException as exc:
                last_err = str(exc)
                time.sleep(3 * (attempt + 1))
        return False, last_err


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def acc_from_id(at_id: str | None) -> str | None:
    if not at_id:
        return None
    parts = [p for p in at_id.split("/") if p]
    return parts[-1] if parts else None


def strip_version(gene_id: str) -> str:
    gid = str(gene_id)
    m = re.match(r"^(ENS[A-Z]*G\d+)(\.\d+)?(_PAR_Y)?$", gid)
    return m.group(1) + (m.group(3) or "") if m else gid


def join(values, sep=";") -> str:
    vals = []
    for v in values:
        if v is None or (isinstance(v, float) and np.isnan(v)):
            continue
        s = str(v)
        if s and s not in vals:
            vals.append(s)
    return sep.join(vals)


def parse_size_range(sr: str | None) -> tuple[float | None, float | None]:
    if not sr:
        return None, None
    s = sr.strip()
    m = re.match(r"^>\s*(\d+)$", s)
    if m:
        return float(m.group(1)), None
    m = re.match(r"^<\s*(\d+)$", s)
    if m:
        return None, float(m.group(1))
    m = re.match(r"^(\d+)\s*-\s*(\d+)$", s)
    if m:
        return float(m.group(1)), float(m.group(2))
    return None, None


def summarize_audit(audit: dict | None) -> tuple[str, str, int, int, int]:
    if not audit:
        return "", "", 0, 0, 0
    parts = []
    for level in ["ERROR", "NOT_COMPLIANT", "WARNING", "INTERNAL_ACTION"]:
        cats = sorted({a.get("category", "") for a in audit.get(level, [])})
        if cats:
            parts.append(f"{level}: " + ", ".join(cats))
    worst = next((lv for lv in ["ERROR", "NOT_COMPLIANT", "WARNING"] if audit.get(lv)), "")
    return (" | ".join(parts), worst, len(audit.get("ERROR", [])),
            len(audit.get("NOT_COMPLIANT", [])), len(audit.get("WARNING", [])))


# ---------------------------------------------------------------------------
# Metadata retrieval
# ---------------------------------------------------------------------------


def fetch_metadata(client: EncodeClient) -> dict:
    md: dict = {}
    assay_params = [("assay_title", a) for a in ASSAY_TITLES_CONSIDERED]

    log.info("Querying experiments ...")
    md["experiments"] = client.search(
        [("type", "Experiment"), ("status", "released"),
         ("replicates.library.biosample.donor.organism.scientific_name", HUMAN)]
        + assay_params + [("frame", "object")],
        "experiments_human_rnaseq_family")
    log.info("  %d experiments", len(md["experiments"]))

    log.info("Querying replicates ...")
    md["replicates"] = client.search(
        [("type", "Replicate")] + [("experiment.assay_title", a) for a in ASSAY_TITLES_CONSIDERED]
        + [("frame", "object")], "replicates_rnaseq_family")
    log.info("  %d replicates", len(md["replicates"]))

    log.info("Querying libraries ...")
    md["libraries"] = client.search(
        [("type", "Library"), ("nucleic_acid_term_name", "RNA"), ("frame", "object")],
        "libraries_rna")
    log.info("  %d libraries", len(md["libraries"]))
    # Libraries referenced by replicates but missing from the RNA-library query (e.g. RNA-seq
    # libraries mis-recorded with nucleic_acid_term_name "DNA") are fetched explicitly.
    have = {l["@id"] for l in md["libraries"]}
    missing = sorted({r.get("library") for r in md["replicates"] if r.get("library") and r.get("library") not in have})
    for i in range(0, len(missing), 80):
        chunk = missing[i:i + 80]
        key = hashlib.sha1(",".join(chunk).encode()).hexdigest()[:16]
        md["libraries"] += client.search([("type", "Library")] + [("accession", acc_from_id(x)) for x in chunk]
                                         + [("frame", "object")], f"libraries_extra_{key}")
    log.info("  %d libraries after fetching %d referenced-but-missing", len(md["libraries"]), len(missing))

    log.info("Querying biosamples ...")
    md["biosamples"] = client.search(
        [("type", "Biosample"), ("organism.scientific_name", HUMAN), ("frame", "object")],
        "biosamples_human")
    log.info("  %d biosamples", len(md["biosamples"]))

    md["biosample_types"] = client.search([("type", "BiosampleType"), ("frame", "object")], "biosample_types")
    md["donors"] = client.search([("type", "HumanDonor"), ("frame", "object")], "human_donors")
    md["treatments"] = client.search([("type", "Treatment"), ("frame", "object")], "treatments")
    md["platforms"] = client.search([("type", "Platform"), ("frame", "object")], "platforms")
    md["labs"] = client.search([("type", "Lab"), ("frame", "object")], "labs")

    log.info("Querying gene-quantification files ...")
    md["gene_quant_files"] = client.search(
        [("type", "File"), ("output_type", "gene quantifications")] + assay_params
        + [("frame", "object")], "files_gene_quantifications")
    log.info("  %d gene quantification files (all statuses)", len(md["gene_quant_files"]))

    log.info("Querying fastq files ...")
    fastq_fields = ["accession", "dataset", "biological_replicates", "technical_replicates",
                    "replicate", "read_length", "run_type", "paired_end", "platform", "status",
                    "read_count", "mapped_read_length"]
    md["fastq_files"] = client.search(
        [("type", "File"), ("file_format", "fastq"), ("assay_title", PRIMARY_TRACK_ASSAY),
         ("assay_title", SMALL_TRACK_ASSAY)] + [("field", f) for f in fastq_fields],
        "files_fastq_total_small")
    log.info("  %d fastq files", len(md["fastq_files"]))

    log.info("Querying reference GTFs ...")
    md["reference_gtfs"] = client.search(
        [("type", "File"), ("file_format", "gtf"), ("output_type", "genome reference"),
         ("frame", "object")], "reference_gtfs")
    return md


def fetch_analyses(client: EncodeClient, accessions: list[str]) -> dict[str, dict]:
    out: dict[str, dict] = {}
    accs = sorted(set(a for a in accessions if a))
    batch = 80
    for i in range(0, len(accs), batch):
        chunk = accs[i:i + batch]
        key = hashlib.sha1(",".join(chunk).encode()).hexdigest()[:16]
        res = client.search([("type", "Analysis")] + [("accession", a) for a in chunk]
                            + [("frame", "object")], f"analyses_{key}")
        for a in res:
            out[a["accession"]] = a
    return out


def fetch_documents(client: EncodeClient, doc_ids: list[str]) -> dict[str, dict]:
    out: dict[str, dict] = {}
    ids = sorted(set(doc_ids))
    batch = 60
    for i in range(0, len(ids), batch):
        chunk = ids[i:i + batch]
        key = hashlib.sha1(",".join(chunk).encode()).hexdigest()[:16]
        res = client.search([("type", "Document")] + [("uuid", acc_from_id(d)) for d in chunk]
                            + [("frame", "object")], f"documents_{key}")
        for d in res:
            out[d["@id"]] = d
    return out


# ---------------------------------------------------------------------------
# Classification
# ---------------------------------------------------------------------------


def library_rna_class(lib: dict, assay_title: str) -> tuple[str, str]:
    """Return (library_rna_class, vault_rna_suitability) from library metadata."""
    dep = set(lib.get("depleted_in_term_name") or [])
    lo, hi = parse_size_range(lib.get("size_range"))
    sr = lib.get("size_range") or "not recorded"
    if "polyadenylated mRNA" in dep:
        cls = "polyA-depleted RNA"
    elif assay_title == "polyA plus RNA-seq":
        cls = "polyA-selected RNA"
    elif hi is not None and hi <= 200 and lo is None:
        cls = f"small RNA fraction (size_range {sr})"
    elif "rRNA" in dep:
        cls = f"rRNA-depleted total RNA (size_range {sr})"
    else:
        cls = f"total RNA, depletion not recorded (size_range {sr})"

    if assay_title == SMALL_TRACK_ASSAY or (hi is not None and hi <= 200 and lo is None):
        suit = "most suitable: <200 nt fraction retains ~100 nt vault RNAs"
    elif lo is not None and lo >= 200 and hi is None:
        suit = "poor: RNA size-selected >200 nt; ~100 nt vault RNAs largely excluded"
    elif lo is not None and hi is not None:
        suit = (f"poor/uncertain: insert/fragment size range {sr} recorded; mature vault RNAs "
                "(88-100 nt) shorter than typical inserts")
    else:
        suit = "uncertain: no RNA size selection recorded"
    return cls, suit


def classify_experiment(exp: dict, ctx: dict) -> dict:
    """Decide track membership and eligibility for one experiment."""
    reasons: list[str] = []
    flags: list[str] = []
    assay_title = exp.get("assay_title")
    bt = ctx["biosample_types"].get(exp.get("biosample_ontology"), {})
    classification = bt.get("classification", "")

    if assay_title == PRIMARY_TRACK_ASSAY:
        track = "total_RNA"
    elif assay_title == SMALL_TRACK_ASSAY:
        track = "small_RNA"
    else:
        track = "none"
        reasons.append(f"assay_title '{assay_title}': " + ASSAY_EXCLUSION_REASON.get(
            assay_title, "not total RNA-seq (or the separately reported small RNA-seq track)"))

    if classification not in ELIGIBLE_CLASSIFICATIONS:
        reasons.append(f"biosample classification '{classification}' is not primary cell/tissue")

    if exp.get("perturbed"):
        reasons.append("experiment flagged perturbed=true by ENCODE")

    reps = ctx["reps_by_exp"].get(exp["@id"], [])
    if not reps:
        reasons.append("no released replicates")

    biosamples, libraries = [], []
    for rep in reps:
        lib = ctx["libraries"].get(rep.get("library"))
        if not lib:
            reasons.append(f"library metadata for replicate {rep.get('biological_replicate_number')}_{rep.get('technical_replicate_number')} could not be retrieved")
            continue
        libraries.append(lib)
        if lib.get("nucleic_acid_term_name") != "RNA":
            flags.append(f"library {lib['accession']} nucleic_acid_term_name recorded as '{lib.get('nucleic_acid_term_name')}' (ENCODE metadata inconsistency for an RNA-seq library)")
        bs = ctx["biosamples"].get(lib.get("biosample"))
        if bs:
            biosamples.append(bs)
        else:
            reasons.append(f"biosample metadata for library {lib['accession']} could not be retrieved")

    for bs in biosamples:
        bst = ctx["biosample_types"].get(bs.get("biosample_ontology"), {})
        if bst.get("classification") and bst.get("classification") not in ELIGIBLE_CLASSIFICATIONS:
            reasons.append(f"replicate biosample {bs['accession']} classification '{bst.get('classification')}'")
        dis = " ".join(bs.get("disease_term_name") or []) if isinstance(bs.get("disease_term_name"), list) else str(bs.get("disease_term_name") or "")
        text = " ".join([dis, bst.get("term_name", ""), bs.get("summary", ""), bs.get("description", "") or ""])
        if CANCER_PATTERN.search(dis) or CANCER_PATTERN.search(bst.get("term_name", "")):
            reasons.append(f"tumour/cancer biosample {bs['accession']} ({dis or bst.get('term_name')})")
        if IMMORTAL_PATTERN.search(text):
            reasons.append(f"biosample {bs['accession']} described as immortalized/transformed/organoid/xenograft/iPSC")
        if IN_VITRO_STATE_PATTERN.search(bst.get("term_name", "")):
            reasons.append(f"biosample term '{bst.get('term_name')}' denotes an in vitro activated/stimulated/polarized state "
                           "(no untreated baseline; excluded as experimentally manipulated even where ENCODE lists no treatment)")
        elif POSSIBLE_POLARIZED_PATTERN.search(bst.get("term_name", "")):
            flags.append(f"biosample term '{bst.get('term_name')}' may reflect ex vivo sorting or in vitro polarization; ENCODE records no treatment")
        if bs.get("host"):
            reasons.append(f"biosample {bs['accession']} has a host organism (xenograft)")
        if bs.get("treatments"):
            tnames = [ctx["treatments"].get(t, {}).get("treatment_term_name", acc_from_id(t)) for t in bs["treatments"]]
            reasons.append(f"biosample {bs['accession']} treated in vitro ({join(tnames, ', ')})")
        if bs.get("genetic_modifications") or bs.get("applied_modifications"):
            reasons.append(f"biosample {bs['accession']} genetically modified")
        if bs.get("perturbed"):
            reasons.append(f"biosample {bs['accession']} flagged perturbed")
        if dis:
            flags.append(f"donor/biosample disease: {dis}")
    for lib in libraries:
        if lib.get("treatments"):
            # Library treatments in ENCODE are enzymatic library-prep steps (e.g. Tobacco Acid
            # Pyrophosphatase, Antarctic Phosphatase, DSN normalization), not biological
            # perturbations: record, do not exclude.
            tnames = [ctx["treatments"].get(t, {}).get("treatment_term_name", acc_from_id(t)) for t in lib["treatments"]]
            flags.append(f"library prep enzymatic treatment: {join(tnames, ', ')}")
            if any("DSN" in (t or "") or "duplex-specific nuclease" in (t or "") for t in tnames):
                flags.append("DSN normalization applied: relative abundances distorted")
        if track == "total_RNA" and "polyadenylated mRNA" in (lib.get("depleted_in_term_name") or []):
            reasons.append(f"library {lib['accession']} depleted of polyA RNA (not total RNA)")
        lo, hi = parse_size_range(lib.get("size_range"))
        if track == "total_RNA" and hi is not None and lo is None and hi <= 200:
            flags.append(f"library {lib['accession']} size_range {lib.get('size_range')} looks like a small-RNA fraction")

    # deduplicate reasons
    seen, uniq = set(), []
    for r in reasons:
        if r not in seen:
            uniq.append(r)
            seen.add(r)
    return {"track": track, "eligible": track != "none" and not uniq,
            "reasons": uniq, "flags": sorted(set(flags))}


# ---------------------------------------------------------------------------
# Quantification-file selection
# ---------------------------------------------------------------------------

ANNOT_RANK = {"V19": 19, "V22": 22, "V24": 24, "V29": 29, "V40": 40}


def select_quant_files(exp: dict, files: list[dict], analyses: dict) -> tuple[list[dict], list[tuple[dict, str]], str]:
    """Select one released gene-quantification TSV per biological replicate.

    Returns (selected_files, [(not_selected_file, reason)], selection_note).
    """
    released = [f for f in files if f.get("status") == "released" and f.get("file_format") == "tsv"]
    not_sel = [(f, f"file status '{f.get('status')}' (only released files used)") for f in files
               if f.get("status") != "released"]
    not_sel += [(f, f"file_format '{f.get('file_format')}'") for f in files
                if f.get("status") == "released" and f.get("file_format") != "tsv"]
    if not released:
        return [], not_sel, "no released gene quantification TSV"

    def group_key(f):
        an = f.get("analyses") or []
        return acc_from_id(an[0]) if an else f"noanalysis:{f.get('assembly')}:{f.get('genome_annotation')}:{acc_from_id(f.get('lab'))}"

    groups = defaultdict(list)
    for f in released:
        groups[group_key(f)].append(f)

    default = acc_from_id(exp.get("default_analysis"))
    note = ""
    if default in groups:
        chosen = default
        note = f"default analysis {default}"
    else:
        def score(k):
            fs = groups[k]
            f0 = fs[0]
            an = analyses.get(k, {})
            uniform = 1 if an.get("pipeline_award_rfas") or (an.get("lab") or "").endswith("encode-processing-pipeline/") else 0
            return (1 if f0.get("assembly") == "GRCh38" else 0,
                    ANNOT_RANK.get(f0.get("genome_annotation") or "", 0),
                    uniform, max(x.get("date_created", "") for x in fs))
        chosen = max(groups, key=score)
        note = (f"default analysis {default} has no released gene quantifications; " if default else
                "no default analysis; ") + f"chose {chosen} (GRCh38 > newest GENCODE > uniform pipeline > newest)"

    for k, fs in groups.items():
        if k != chosen:
            f0 = fs[0]
            for f in fs:
                not_sel.append((f, f"non-selected analysis {k} ({f0.get('assembly')} {f0.get('genome_annotation')}); selected {chosen}"))

    by_rep = defaultdict(list)
    for f in groups[chosen]:
        br = tuple(sorted(f.get("biological_replicates") or []))
        by_rep[br].append(f)

    selected = []
    for br, fs in sorted(by_rep.items()):
        if len(br) != 1:
            for f in fs:
                not_sel.append((f, f"covers biological replicates {list(br)} (pooled), not a single biological replicate"))
            continue
        fs_sorted = sorted(fs, key=lambda f: (bool(f.get("preferred_default")), f.get("date_created", "")), reverse=True)
        selected.append(fs_sorted[0])
        for f in fs_sorted[1:]:
            not_sel.append((f, f"duplicate output for biological replicate {br[0]} in analysis {chosen}; selected {fs_sorted[0]['accession']} (preferred_default/newest)"))
    return selected, not_sel, note


# ---------------------------------------------------------------------------
# GTF gene-ID resolution
# ---------------------------------------------------------------------------


def choose_reference_gtf(ref_gtfs: list[dict], assembly: str, annotation: str) -> dict | None:
    cands = [f for f in ref_gtfs if f.get("assembly") == assembly and f.get("genome_annotation") == annotation]
    if not cands:
        return None

    def score(f):
        name = (f.get("submitted_file_name") or "").lower()
        return (f.get("status") == "released",
                "trna" in name and "ercc" in name,          # merged annotation used by RNA-seq pipelines
                "mirna" not in name and "promoter" not in name and "trnas.gtf" not in name,
                f.get("file_size") or 0)
    best = max(cands, key=score)
    name = (best.get("submitted_file_name") or "").lower()
    if "mirna" in name or "promoter" in name:
        return None
    return best


def resolve_gene_ids(client: EncodeClient, paths: Paths, ref: dict, label: str) -> list[dict]:
    """Return gene records in the GTF whose gene_name is VTRNA1-* or id matches a seed."""
    acc = ref.get("accession") or acc_from_id(ref["@id"])
    dest = os.path.join(paths.refs, f"{acc}.gtf.gz")
    if not os.path.exists(dest):
        ok, msg = client.download(ref["href"], dest, ref.get("md5sum"))
        if not ok:
            raise RuntimeError(f"Failed to download reference GTF {acc}: {msg}")
    seeds = set(SEED_STABLE_IDS.values())
    out = []
    name_re = re.compile(r'gene_name "([^"]+)"')
    id_re = re.compile(r'gene_id "([^"]+)"')
    type_re = re.compile(r'gene_(?:bio)?type "([^"]+)"')
    seen = set()
    with gzip.open(dest, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            if "VTRNA1" not in line and not any(s in line for s in seeds):
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 9:
                continue
            attrs = cols[8]
            gname = (name_re.search(attrs) or [None, None])[1]
            gid = (id_re.search(attrs) or [None, None])[1]
            if gid is None:
                continue
            if not ((gname and gname.startswith("VTRNA1-")) or strip_version(gid) in seeds):
                continue
            if gid in seen:
                continue
            if cols[2] not in ("gene", "transcript", "exon"):
                continue
            seen.add(gid)
            out.append({"annotation_label": label, "gtf_accession": acc,
                        "gtf_file_name": ref.get("submitted_file_name"),
                        "gene_symbol": gname, "gene_id_gtf": gid, "stable_gene_id": strip_version(gid),
                        "gene_type": (type_re.search(attrs) or [None, None])[1],
                        "chrom": cols[0], "start": int(cols[3]), "end": int(cols[4]), "strand": cols[6],
                        "gtf_feature_used": cols[2]})
    return out


# ---------------------------------------------------------------------------
# Quantification parsing
# ---------------------------------------------------------------------------


def parse_quant(path_gz: str) -> tuple[str, pd.DataFrame, dict]:
    """Parse an ENCODE gene-quantification TSV.

    Returns (schema, dataframe, file_level_stats).  dataframe always has a
    'gene_id' column (original identifier) plus schema-specific columns.
    """
    with gzip.open(path_gz, "rt") as fh:
        text = fh.read()
    first = text.split("\n", 1)[0]
    stats: dict = {}
    if first.startswith("N_unmapped"):
        df = pd.read_csv(io.StringIO(text), sep="\t", header=None,
                         names=["gene_id", "count_unstranded", "count_strand1", "count_strand2"],
                         dtype={"gene_id": str})
        special = df["gene_id"].str.startswith("N_")
        for _, r in df[special].iterrows():
            stats[r["gene_id"]] = r["count_unstranded"]
        genes = df[~special].copy()
        s1, s2 = genes["count_strand1"].sum(), genes["count_strand2"].sum()
        if s1 > 3 * s2:
            col = "count_strand1"
        elif s2 > 3 * s1:
            col = "count_strand2"
        else:
            col = "count_unstranded"
        genes["count_sense"] = genes[col]
        tot = genes[col].sum()
        genes["CPM"] = genes[col] / tot * 1e6 if tot > 0 else np.nan
        stats.update({"sense_column": col, "total_assigned_counts": int(tot),
                      "sum_strand1": int(s1), "sum_strand2": int(s2), "n_genes": len(genes)})
        return "STAR_ReadsPerGene", genes, stats

    df = pd.read_csv(io.StringIO(text), sep="\t", dtype={"gene_id": str, "transcript_id(s)": str},
                     low_memory=False)
    cols = set(df.columns)
    if {"gene_id", "expected_count", "TPM"} <= cols:
        schema = "RSEM_genes_results"
    else:
        idcol = next((c for c in ["gene_id", "tracking_id", "Geneid", "gene"] if c in cols), None)
        if idcol is None:
            raise ValueError(f"Unrecognised quantification schema: {list(df.columns)[:10]}")
        df = df.rename(columns={idcol: "gene_id"})
        schema = "other:" + ",".join(list(df.columns)[:8])
    for c in ["TPM", "FPKM", "expected_count", "length", "effective_length"]:
        if c in df.columns:
            df[c] = pd.to_numeric(df[c], errors="coerce")
    stats["n_genes"] = len(df)
    if "TPM" in df.columns:
        stats["sum_TPM"] = float(df["TPM"].sum())
    if "expected_count" in df.columns:
        stats["sum_expected_count"] = float(df["expected_count"].sum())
    if "FPKM" in df.columns:
        stats["sum_FPKM"] = float(df["FPKM"].sum())
    return schema, df, stats


def extract_targets(schema: str, df: pd.DataFrame, gene_map: dict[str, list[dict]],
                    all_ids_by_symbol: dict[str, set]) -> dict[str, dict]:
    """Find VTRNA1 rows in a parsed quantification table."""
    stripped = df["gene_id"].astype(str).map(strip_version)
    symcol = next((c for c in ["gene_name", "gene_short_name", "gene_symbol"] if c in df.columns), None)
    res = {}
    for sym in TARGET_SYMBOLS:
        recs = gene_map.get(sym, [])
        ids_this_annot = {r["stable_gene_id"] for r in recs}
        ids_any = set(all_ids_by_symbol.get(sym, set())) | ids_this_annot
        if sym in SEED_STABLE_IDS:
            ids_any.add(SEED_STABLE_IDS[sym])
        mask = stripped.isin(ids_any)
        match_by = "ensembl_id(version-stripped)"
        if symcol is not None:
            smask = df[symcol].astype(str) == sym
            if smask.any() and not mask.any():
                match_by = "gene_symbol_column"
            mask = mask | smask
        hits = df[mask]
        out = {"gtf_gene_ids": join(r["gene_id_gtf"] for r in recs),
               "symbol_in_annotation": bool(recs)}
        if hits.empty:
            out["status"] = "absent_from_file" if recs else "absent_from_annotation_and_file"
            res[sym] = out
            continue
        if len(hits) > 1:
            # Prefer exact ID from this annotation; keep PAR copies out
            prim = hits[~hits["gene_id"].astype(str).str.contains("_PAR_Y")]
            hits = prim if not prim.empty else hits
        row = hits.iloc[0]
        out.update({"status": "present", "n_matching_rows": int(len(hits)),
                    "gene_id_reported": row["gene_id"], "matched_by": match_by,
                    "id_in_this_annotation": strip_version(row["gene_id"]) in ids_this_annot,
                    "gene_symbol_reported_in_file": row[symcol] if symcol else None})
        if schema == "STAR_ReadsPerGene":
            for c in ["count_unstranded", "count_strand1", "count_strand2", "count_sense", "CPM"]:
                out[c] = float(row[c])
        else:
            for c in ["expected_count", "TPM", "FPKM", "length", "effective_length",
                      "posterior_mean_count", "pme_TPM", "pme_FPKM"]:
                if c in df.columns:
                    out[c] = float(pd.to_numeric(row[c], errors="coerce"))
            if "transcript_id(s)" in df.columns:
                out["transcript_ids"] = row["transcript_id(s)"]
        res[sym] = out
    return res


# ---------------------------------------------------------------------------
# Main pipeline
# ---------------------------------------------------------------------------


def build_context(md: dict) -> dict:
    bts = {b["@id"]: b for b in md["biosample_types"]}
    # Some BiosampleType objects are not returned by the bulk search (e.g. newer NTR
    # terms); derive the classification from the @id name, which encodes it.
    referenced = {e.get("biosample_ontology") for e in md["experiments"]} | \
                 {b.get("biosample_ontology") for b in md["biosamples"]}
    for ref in referenced:
        if ref and ref not in bts:
            name = acc_from_id(ref)
            for cls in ["in_vitro_differentiated_cells", "primary_cell", "cell_line", "tissue",
                        "organoid", "whole_organisms", "cell_free_sample"]:
                if name.startswith(cls + "_"):
                    bts[ref] = {"@id": ref, "classification": cls.replace("_", " ").replace("in vitro", "in vitro"),
                                "term_id": name[len(cls) + 1:].replace("_", ":", 1), "term_name": "",
                                "_derived_from_id": True}
                    break
    ctx = {
        "biosample_types": bts,
        "libraries": {l["@id"]: l for l in md["libraries"]},
        "biosamples": {b["@id"]: b for b in md["biosamples"]},
        "donors": {d["@id"]: d for d in md["donors"]},
        "treatments": {t["@id"]: t for t in md["treatments"]},
        "platforms": {p["@id"]: p for p in md["platforms"]},
        "labs": {l["@id"]: l for l in md["labs"]},
    }
    reps_by_exp = defaultdict(list)
    for r in md["replicates"]:
        if r.get("status") in ("released",):
            reps_by_exp[r["experiment"]].append(r)
    ctx["reps_by_exp"] = reps_by_exp
    files_by_exp = defaultdict(list)
    for f in md["gene_quant_files"]:
        f["_annotation_source"] = "file.genome_annotation" if f.get("genome_annotation") else ""
        if not f.get("genome_annotation"):
            mm = re.search(r"GENCODE\s*v?(\d+)", f.get("submitter_comment") or "", re.I)
            if mm:
                f["genome_annotation"] = f"V{mm.group(1)}"
                f["_annotation_source"] = "file.submitter_comment"
            else:
                f["_annotation_source"] = "not recorded"
        files_by_exp[f["dataset"]].append(f)
    ctx["gq_by_exp"] = files_by_exp
    fq = defaultdict(list)
    for f in md["fastq_files"]:
        if f.get("status") in ("released",):
            ds = f["dataset"]["@id"] if isinstance(f.get("dataset"), dict) else f.get("dataset")
            for br in f.get("biological_replicates") or []:
                fq[(ds, br)].append(f)
    ctx["fastq_by_exp_rep"] = fq
    return ctx


def replicate_metadata(exp: dict, br: int, ctx: dict, documents: dict, file_techreps: list | None = None) -> dict:
    all_reps = [r for r in ctx["reps_by_exp"].get(exp["@id"], []) if r.get("biological_replicate_number") == br]
    # restrict to the technical replicate(s) actually contained in the selected file
    tr = set(file_techreps or [])
    reps = [r for r in all_reps if f"{br}_{r.get('technical_replicate_number')}" in tr] or all_reps
    libs = [ctx["libraries"].get(r.get("library"), {}) for r in reps]
    bss = [ctx["biosamples"].get(l.get("biosample"), {}) for l in libs if l]
    bts = [ctx["biosample_types"].get(b.get("biosample_ontology"), {}) for b in bss if b]
    donors = [ctx["donors"].get(b.get("donor"), {}) for b in bss if b]
    fqs = ctx["fastq_by_exp_rep"].get((exp["@id"], br), [])
    if tr:
        fq_tr = [f for f in fqs if set(f.get("technical_replicates") or []) & tr]
        fqs = fq_tr or fqs
    classes, suits = zip(*[library_rna_class(l, exp.get("assay_title")) for l in libs if l]) if any(libs) else ((), ())
    docs = []
    for l in libs:
        for d in l.get("documents") or []:
            doc = documents.get(d, {})
            att = (doc.get("attachment") or {}).get("download") or ""
            docs.append(f"{doc.get('document_type', '')}:{att}".strip(":"))
    platforms = []
    for f in fqs:
        pl = f.get("platform")
        if isinstance(pl, dict):
            platforms.append(pl.get("term_name") or acc_from_id(pl.get("@id")))
        elif pl:
            platforms.append(ctx["platforms"].get(pl, {}).get("term_name", acc_from_id(pl)))
    bt0 = bts[0] if bts else {}
    disease = []
    for b in bss:
        d = b.get("disease_term_name")
        disease += d if isinstance(d, list) else ([d] if d else [])
    return {
        "biosample_accession": join(b.get("accession") for b in bss),
        "biosample_term_name": join(bt.get("term_name") for bt in bts),
        "biosample_term_id": join(bt.get("term_id") for bt in bts),
        "biosample_type": join(bt.get("classification") for bt in bts),
        "organ": join(o for bt in bts for o in (bt.get("organ_slims") or [])),
        "system": join(o for bt in bts for o in (bt.get("system_slims") or [])),
        "cell_type": join(o for bt in bts for o in (bt.get("cell_slims") or [])) if bt0.get("classification") == "primary cell" else "",
        "cell_slims": join(o for bt in bts for o in (bt.get("cell_slims") or [])),
        "anatomical_site": join(bt.get("term_name") for bt in bts) if bt0.get("classification") == "tissue" else "",
        "biosample_summary": join((b.get("summary") for b in bss), " || "),
        "biosample_description": join((b.get("description") for b in bss), " || "),
        "preservation_method": join(b.get("preservation_method") for b in bss),
        "health_status": join((b.get("health_status") for b in bss), " || "),
        "disease_term_name": join(disease),
        "donor_accession": join(d.get("accession") for d in donors),
        "donor_sex": join(b.get("sex") for b in bss) or join(d.get("sex") for d in donors),
        "donor_age": join(b.get("age") for b in bss) or join(d.get("age") for d in donors),
        "donor_age_units": join(b.get("age_units") for b in bss) or join(d.get("age_units") for d in donors),
        "donor_life_stage": join(b.get("life_stage") for b in bss) or join(d.get("life_stage") for d in donors),
        "donor_ethnicity": join(e for d in donors for e in (d.get("ethnicity") or [])),
        "technical_replicates": join(f"{br}_{r.get('technical_replicate_number')}" for r in reps),
        "n_technical_replicates": len(reps),
        "all_technical_replicates_of_bioreplicate": join(f"{br}_{r.get('technical_replicate_number')}" for r in all_reps),
        "library_accession": join(l.get("accession") for l in libs),
        "library_rna_class": join(classes),
        "vault_rna_suitability": join(suits),
        "library_size_range": join(l.get("size_range") for l in libs),
        "library_depleted_in": join(x for l in libs for x in (l.get("depleted_in_term_name") or [])),
        "library_size_selection_method": join(l.get("library_size_selection_method") for l in libs),
        "library_extraction_method": join(l.get("extraction_method") for l in libs),
        "library_fragmentation_methods": join(x for l in libs for x in (l.get("fragmentation_methods") or [])),
        "library_average_fragment_size": join(l.get("average_fragment_size") for l in libs),
        "library_spikeins": join(acc_from_id(x) for l in libs for x in (l.get("spikeins_used") or [])),
        "library_treatments": join(ctx["treatments"].get(t, {}).get("treatment_term_name", acc_from_id(t))
                                   for l in libs for t in (l.get("treatments") or [])),
        "library_rna_integrity_number": join(l.get("rna_integrity_number") for l in libs),
        "library_documents": join(docs),
        "strandedness": join(l.get("strand_specificity") for l in libs),
        "read_length": join(sorted({f.get("read_length") for f in fqs if f.get("read_length")})),
        "run_type": join(f.get("run_type") for f in fqs),
        "sequencing_platform": join(platforms),
        "n_fastq_files": len(fqs),
    }


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--outdir", default=os.path.dirname(os.path.abspath(__file__)))
    ap.add_argument("--refresh", action="store_true", help="re-query the ENCODE API instead of using cached responses")
    ap.add_argument("--workers", type=int, default=6)
    ap.add_argument("--limit", type=int, default=0, help="debug: process only the first N eligible experiments")
    args = ap.parse_args(argv)

    paths = Paths(args.outdir)
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s",
                        handlers=[logging.StreamHandler(sys.stderr),
                                  logging.FileHandler(os.path.join(paths.logs, "run.log"))])
    run_started = utcnow()
    client = EncodeClient(paths, refresh=args.refresh)

    md = fetch_metadata(client)
    ctx = build_context(md)

    experiments = sorted(md["experiments"], key=lambda e: e["accession"])
    # experiment-level classification
    classified = {e["@id"]: classify_experiment(e, ctx) for e in experiments}

    eligible = [e for e in experiments if classified[e["@id"]]["eligible"]]
    if args.limit:
        eligible = eligible[: args.limit]
    log.info("Eligible experiments: %d (total_RNA %d, small_RNA %d)", len(eligible),
             sum(classified[e["@id"]]["track"] == "total_RNA" for e in eligible),
             sum(classified[e["@id"]]["track"] == "small_RNA" for e in eligible))

    # analyses for eligible experiments (default + all referenced by quant files)
    an_accs = []
    for e in eligible:
        an_accs += [acc_from_id(a) for a in e.get("analyses") or []]
    analyses = fetch_analyses(client, an_accs)
    log.info("Fetched %d analyses", len(analyses))

    doc_ids = []
    for e in eligible:
        for r in ctx["reps_by_exp"].get(e["@id"], []):
            doc_ids += ctx["libraries"].get(r.get("library"), {}).get("documents") or []
    documents = fetch_documents(client, doc_ids)

    # select files
    selections = []      # (exp, file, track, note)
    manifest_rows = []
    for e in experiments:
        c = classified[e["@id"]]
        bt = ctx["biosample_types"].get(e.get("biosample_ontology"), {})
        base = {
            "experiment_accession": e["accession"], "assay_title": e.get("assay_title"),
            "track": c["track"], "biosample_term_name": bt.get("term_name") or e.get("biosample_summary"),
            "biosample_type": bt.get("classification"), "lab": ctx["labs"].get(e.get("lab"), {}).get("title", acc_from_id(e.get("lab"))),
            "perturbed": e.get("perturbed"), "experiment_status": e.get("status"),
            "default_analysis": acc_from_id(e.get("default_analysis")),
            "experiment_flags": "; ".join(c["flags"]),
        }
        audit_txt, worst, n_err, n_nc, n_warn = summarize_audit(e.get("audit"))
        base.update({"experiment_audit_worst": worst, "experiment_audit": audit_txt})
        if not c["eligible"]:
            manifest_rows.append({**base, "row_type": "experiment_excluded", "included": False,
                                  "exclusion_reason": "; ".join(c["reasons"])})
            continue
        if e not in eligible:   # --limit debug mode
            continue
        files = ctx["gq_by_exp"].get(e["@id"], [])
        sel, not_sel, note = select_quant_files(e, files, analyses)
        if not sel:
            manifest_rows.append({**base, "row_type": "eligible_no_gene_quantification", "included": False,
                                  "exclusion_reason": "experiment meets biosample/assay criteria but has no released single-replicate gene quantification TSV"
                                  + (f" ({note})" if note else "")})
        for f in sel:
            selections.append((e, f, c["track"], note))
            an = analyses.get(acc_from_id((f.get("analyses") or [None])[0]), {})
            manifest_rows.append({**base, "row_type": "selected_file", "included": True, "exclusion_reason": "",
                                  "file_accession": f["accession"], "biological_replicate": join(f.get("biological_replicates")),
                                  "technical_replicates": join(f.get("technical_replicates")),
                                  "file_status": f.get("status"), "genome_assembly": f.get("assembly"),
                                  "annotation_version": f.get("genome_annotation"),
                                  "annotation_version_source": f.get("_annotation_source"),
                                  "analysis_accession": an.get("accession"), "analysis_title": an.get("title"),
                                  "selection_note": note})
        for f, why in not_sel:
            an = analyses.get(acc_from_id((f.get("analyses") or [None])[0]), {})
            manifest_rows.append({**base, "row_type": "alternate_file_not_selected", "included": False,
                                  "exclusion_reason": why, "file_accession": f["accession"],
                                  "biological_replicate": join(f.get("biological_replicates")),
                                  "technical_replicates": join(f.get("technical_replicates")),
                                  "file_status": f.get("status"), "genome_assembly": f.get("assembly"),
                                  "annotation_version": f.get("genome_annotation"),
                                  "analysis_accession": an.get("accession") or acc_from_id((f.get("analyses") or [None])[0]),
                                  "analysis_title": an.get("title"), "selection_note": note})

    log.info("Selected %d quantification files", len(selections))

    # resolve gene ids per annotation
    annot_keys = sorted({(f.get("assembly"), f.get("genome_annotation")) for _, f, _, _ in selections}
                        | {("GRCh38", "V29"), ("GRCh38", "V24"), ("hg19", "V19")},
                        key=lambda k: (str(k[0]), str(k[1])))
    gene_maps: dict[tuple, dict[str, list[dict]]] = {}
    gtf_records = []
    for asm, ann in annot_keys:
        ref = choose_reference_gtf(md["reference_gtfs"], asm, ann)
        if ref is None:
            log.warning("No ENCODE reference GTF for %s %s: genes will be matched by version-stripped "
                        "stable Ensembl ID resolved from the other ENCODE GTFs", asm, ann)
            gene_maps[(asm, ann)] = {}
            continue
        recs = resolve_gene_ids(client, paths, ref, f"{asm} {ann}")
        gtf_records += recs
        gm = defaultdict(list)
        for r in recs:
            gm[r["gene_symbol"]].append(r)
        gene_maps[(asm, ann)] = gm
        log.info("GTF %s (%s %s): %s", ref.get("accession") or ref["@id"], asm, ann,
                 {s: [r["gene_id_gtf"] for r in gm.get(s, [])] for s in TARGET_SYMBOLS})
    all_ids_by_symbol = defaultdict(set)
    for r in gtf_records:
        all_ids_by_symbol[r["gene_symbol"]].add(r["stable_gene_id"])
    pd.DataFrame(gtf_records).to_csv(os.path.join(paths.extract, "vtrna1_gtf_records.csv"), index=False)

    # download + parse
    failed = []
    results = {}

    def work(item):
        e, f, track, note = item
        dest = os.path.join(paths.quant, track, f.get("assembly") or "NA", f"{f['accession']}.tsv.gz")
        ok, msg = client.download(f["href"], dest, f.get("md5sum"))
        if not ok:
            return f["accession"], {"download_ok": False, "download_msg": msg}
        try:
            schema, df, stats = parse_quant(dest)
            hits = extract_targets(schema, df, gene_maps.get((f.get("assembly"), f.get("genome_annotation")), {}),
                                   all_ids_by_symbol)
            return f["accession"], {"download_ok": True, "download_msg": msg, "parsed": True, "schema": schema,
                                    "stats": stats, "hits": hits, "local_path": os.path.relpath(dest, paths.out)}
        except Exception as exc:  # noqa: BLE001
            return f["accession"], {"download_ok": True, "download_msg": msg, "parsed": False,
                                    "parse_error": repr(exc), "local_path": os.path.relpath(dest, paths.out)}

    with ThreadPoolExecutor(max_workers=args.workers) as pool:
        futs = [pool.submit(work, it) for it in selections]
        for i, fut in enumerate(as_completed(futs), 1):
            acc, res = fut.result()
            results[acc] = res
            if not res.get("download_ok"):
                failed.append((acc, res.get("download_msg")))
            if i % 50 == 0 or i == len(futs):
                log.info("  processed %d/%d files", i, len(futs))
    with open(os.path.join(paths.logs, "failed_downloads.tsv"), "w") as fh:
        fh.write("file_accession\terror\n")
        for acc, msg in failed:
            fh.write(f"{acc}\t{msg}\n")

    # build tables
    wide_rows, long_rows, raw_rows = [], [], []
    for e, f, track, note in selections:
        br = (f.get("biological_replicates") or [None])[0]
        meta = replicate_metadata(e, br, ctx, documents, f.get("technical_replicates"))
        res = results.get(f["accession"], {})
        an = analyses.get(acc_from_id((f.get("analyses") or [None])[0]), {})
        audit_txt, worst, n_err, n_nc, n_warn = summarize_audit(e.get("audit"))
        qm = an.get("quality_metrics_report") or {}
        depth = [x.get("metric") for x in qm.get("Read depth", []) if br in (x.get("biological_replicates") or [])]
        row = {
            "track": track,
            "biosample_term_name": meta["biosample_term_name"],
            "biosample_type": meta["biosample_type"],
            "organ": meta["organ"], "cell_type": meta["cell_type"],
            "anatomical_site": meta["anatomical_site"], "system": meta["system"],
            "donor_accession": meta["donor_accession"], "donor_sex": meta["donor_sex"],
            "donor_age": meta["donor_age"], "donor_age_units": meta["donor_age_units"],
            "donor_life_stage": meta["donor_life_stage"], "donor_ethnicity": meta["donor_ethnicity"],
            "health_status": meta["health_status"], "disease_term_name": meta["disease_term_name"],
            "experiment_accession": e["accession"], "biological_replicate": br,
            "technical_replicates": meta["technical_replicates"],
            "file_technical_replicates": join(f.get("technical_replicates")),
            "n_technical_replicates": meta["n_technical_replicates"],
            "all_technical_replicates_of_bioreplicate": meta["all_technical_replicates_of_bioreplicate"],
            "other_technical_replicate_files_not_selected": join(
                x["file_accession"] for x in manifest_rows
                if x.get("experiment_accession") == e["accession"] and x.get("row_type") == "alternate_file_not_selected"
                and str(x.get("exclusion_reason", "")).startswith(f"duplicate output for biological replicate {br} ")),
            "biosample_accession": meta["biosample_accession"],
            "library_accession": meta["library_accession"], "file_accession": f["accession"],
            "file_status": f.get("status"), "file_md5sum": f.get("md5sum"),
            "file_date_created": f.get("date_created"),
            "assay_title": e.get("assay_title"),
            "lab": ctx["labs"].get(e.get("lab"), {}).get("title"),
            "library_rna_class": meta["library_rna_class"], "vault_rna_suitability": meta["vault_rna_suitability"],
            "library_size_range": meta["library_size_range"], "library_depleted_in": meta["library_depleted_in"],
            "library_size_selection_method": meta["library_size_selection_method"],
            "library_extraction_method": meta["library_extraction_method"],
            "library_fragmentation_methods": meta["library_fragmentation_methods"],
            "library_average_fragment_size": meta["library_average_fragment_size"],
            "library_rna_integrity_number": meta["library_rna_integrity_number"],
            "library_spikeins": meta["library_spikeins"], "library_treatments": meta["library_treatments"],
            "library_documents": meta["library_documents"],
            "strandedness": meta["strandedness"], "read_length": meta["read_length"],
            "run_type": meta["run_type"], "sequencing_platform": meta["sequencing_platform"],
            "n_fastq_files": meta["n_fastq_files"],
            "genome_assembly": f.get("assembly"), "annotation_version": f.get("genome_annotation"),
            "annotation_version_source": f.get("_annotation_source"),
            "file_submitter_comment": f.get("submitter_comment"),
            "analysis_step_version": acc_from_id(f.get("analysis_step_version")),
            "analysis_accession": an.get("accession"), "analysis_title": an.get("title"),
            "pipeline_version": an.get("pipeline_version"), "is_default_analysis": an.get("accession") == acc_from_id(e.get("default_analysis")),
            "selection_note": note, "read_depth_qc": join(depth),
            "experiment_audit_worst": worst, "experiment_audit_n_error": n_err,
            "experiment_audit_n_not_compliant": n_nc, "experiment_audit_n_warning": n_warn,
            "experiment_audit": audit_txt,
            "experiment_flags": "; ".join(classified[e["@id"]]["flags"]),
            "biosample_summary": meta["biosample_summary"], "preservation_method": meta["preservation_method"],
            "quant_schema": res.get("schema"), "download_ok": res.get("download_ok"),
            "parsed_ok": res.get("parsed", False), "parse_error": res.get("parse_error"),
            "local_file": res.get("local_path"),
        }
        stats = res.get("stats") or {}
        row["file_n_genes"] = stats.get("n_genes")
        row["file_sum_TPM"] = stats.get("sum_TPM")
        row["file_sum_expected_count"] = stats.get("sum_expected_count")
        row["file_total_assigned_counts"] = stats.get("total_assigned_counts")
        row["file_sense_count_column"] = stats.get("sense_column")
        hits = res.get("hits") or {}
        for sym in TARGET_SYMBOLS:
            p = sym.replace("-", "_")
            h = hits.get(sym)
            if not res.get("download_ok"):
                status, reason = "NA_file_unavailable", f"download failed: {res.get('download_msg')}"
            elif not res.get("parsed"):
                status, reason = "NA_parse_failed", res.get("parse_error")
            elif h["status"] != "present":
                status = "NA_" + h["status"]
                reason = ("gene row absent from quantification file although annotation contains "
                          + h.get("gtf_gene_ids", "") if h["status"] == "absent_from_file"
                          else "gene not in annotation nor file")
            else:
                status, reason = "present", ""
            hh = h or {}
            gene_id = hh.get("gene_id_reported") if status == "present" else None
            has_gtf = bool(gene_maps.get((f.get("assembly"), f.get("genome_annotation"))))
            sym_src = ("GENCODE GTF " + (row["genome_assembly"] or "") + " " + (row["annotation_version"] or "")
                       if has_gtf else "no ENCODE GTF for this annotation; matched by stable Ensembl ID resolved from ENCODE GENCODE V19/V24/V29 GTFs")
            rec = {"gene_id": gene_id, "gene_symbol_source": sym_src,
                   "gtf_gene_id": hh.get("gtf_gene_ids"), "measurement_status": status, "na_reason": reason}
            if res.get("schema") != "STAR_ReadsPerGene":
                eff = hh.get("effective_length")
                rec.update({"TPM": hh.get("TPM"), "expected_count": hh.get("expected_count"), "FPKM": hh.get("FPKM"),
                            "length": hh.get("length"), "effective_length": eff,
                            "posterior_mean_count": hh.get("posterior_mean_count"), "pme_TPM": hh.get("pme_TPM"),
                            "count_sense": None, "count_unstranded": None, "count_strand1": None,
                            "count_strand2": None, "CPM": None})
                if status == "present":
                    rec["quantifiable"] = bool(eff is not None and eff > 0)
                    if not rec["quantifiable"]:
                        rec["measurement_status"] = "present_effective_length_0"
                        rec["na_reason"] = ("RSEM effective_length=0 (gene shorter than mean fragment length): "
                                            "TPM/FPKM/expected_count are 0 by construction, not a measured absence")
                else:
                    rec["quantifiable"] = None
            else:
                rec.update({"TPM": None, "expected_count": None, "FPKM": None,
                            "length": None, "effective_length": None, "quantifiable": (True if status == "present" else None),
                            "posterior_mean_count": None, "pme_TPM": None,
                            "count_sense": hh.get("count_sense"), "count_unstranded": hh.get("count_unstranded"),
                            "count_strand1": hh.get("count_strand1"), "count_strand2": hh.get("count_strand2"),
                            "CPM": hh.get("CPM")})
            for k, v in rec.items():
                row[f"{p}_{k}"] = v
            lrow = {"track": track, "quant_schema": res.get("schema"), "biosample": row["biosample_term_name"], "biosample_type": row["biosample_type"],
                    "organ": row["organ"], "cell_type": row["cell_type"],
                    "donor_accession": row["donor_accession"],
                    "experiment_accession": row["experiment_accession"], "biological_replicate": br,
                    "file_accession": f["accession"], "gene_symbol": sym, **rec,
                    "genome_assembly": row["genome_assembly"], "annotation_version": row["annotation_version"],
                    "library_rna_class": row["library_rna_class"], "vault_rna_suitability": row["vault_rna_suitability"]}
            long_rows.append(lrow)
            if status == "present":
                raw_rows.append({"file_accession": f["accession"], "experiment_accession": e["accession"],
                                 "track": track, "gene_symbol": sym, **{k: v for k, v in hh.items()}})
        wide_rows.append(row)

    wide = pd.DataFrame(wide_rows)
    long = pd.DataFrame(long_rows)
    manifest = pd.DataFrame(manifest_rows)
    pd.DataFrame(raw_rows).to_csv(os.path.join(paths.extract, "vtrna1_rows_extracted_from_source_files.csv"), index=False)

    def order_wide(df):
        first = ["biosample_term_name", "biosample_type", "organ", "cell_type", "anatomical_site",
                 "donor_accession", "donor_sex", "donor_age", "donor_age_units", "donor_life_stage",
                 "donor_ethnicity", "experiment_accession", "biological_replicate", "technical_replicates",
                 "biosample_accession", "library_accession", "file_accession", "genome_assembly",
                 "annotation_version"]
        genes = []
        for sym in TARGET_SYMBOLS:
            p = sym.replace("-", "_")
            genes += [c for c in [f"{p}_gene_id", f"{p}_TPM", f"{p}_expected_count", f"{p}_FPKM"] if c in df.columns]
        rest = [c for c in df.columns if c not in first + genes]
        return df[[c for c in first if c in df.columns] + genes + rest]

    if not wide.empty:
        wide["n_rows_same_donor_and_biosample_term"] = wide.groupby(
            ["track", "donor_accession", "biosample_term_name"], dropna=False)["file_accession"].transform("count")
        wide["n_rows_same_donor_and_biosample_term_this_schema"] = wide.groupby(
            ["track", "quant_schema", "donor_accession", "biosample_term_name"], dropna=False)["file_accession"].transform("count")
    tot = order_wide(wide[wide["track"] == "total_RNA"].copy()) if not wide.empty else wide
    sml = order_wide(wide[wide["track"] == "small_RNA"].copy()) if not wide.empty else wide
    tot = tot.sort_values(["biosample_type", "biosample_term_name", "experiment_accession", "biological_replicate"])
    sml = sml.sort_values(["biosample_type", "biosample_term_name", "experiment_accession", "biological_replicate"])
    tot.to_csv(os.path.join(paths.out, "encode_primary_human_vtrna1_expression.csv"), index=False)
    sml.to_csv(os.path.join(paths.out, "encode_primary_human_vtrna1_smallRNA_expression.csv"), index=False)

    long_cols_total = ["biosample", "biosample_type", "organ", "cell_type", "donor_accession",
                       "experiment_accession", "biological_replicate", "file_accession", "gene_symbol",
                       "gene_id", "TPM", "expected_count", "FPKM", "length", "effective_length", "quantifiable",
                       "quant_schema", "count_sense", "CPM", "measurement_status", "na_reason", "gtf_gene_id", "gene_symbol_source",
                       "genome_assembly", "annotation_version", "library_rna_class", "vault_rna_suitability"]
    long_cols_small = ["biosample", "biosample_type", "organ", "cell_type", "donor_accession",
                       "experiment_accession", "biological_replicate", "file_accession", "gene_symbol",
                       "gene_id", "TPM", "expected_count", "FPKM", "quant_schema", "count_sense", "CPM", "count_unstranded",
                       "count_strand1", "count_strand2", "measurement_status", "na_reason", "gtf_gene_id",
                       "gene_symbol_source", "genome_assembly", "annotation_version", "library_rna_class",
                       "vault_rna_suitability"]
    lt = long[long["track"] == "total_RNA"]
    ls = long[long["track"] == "small_RNA"]
    lt[[c for c in long_cols_total if c in lt.columns]].to_csv(
        os.path.join(paths.out, "encode_primary_human_vtrna1_expression_long.csv"), index=False)
    ls[[c for c in long_cols_small if c in ls.columns]].to_csv(
        os.path.join(paths.out, "encode_primary_human_vtrna1_smallRNA_expression_long.csv"), index=False)

    man_cols = ["experiment_accession", "assay_title", "track", "row_type", "included", "exclusion_reason",
                "file_accession", "biological_replicate", "technical_replicates", "file_status",
                "genome_assembly", "annotation_version", "annotation_version_source", "analysis_accession",
                "analysis_title", "default_analysis", "selection_note", "biosample_term_name", "biosample_type", "lab",
                "perturbed", "experiment_status", "experiment_flags", "experiment_audit_worst", "experiment_audit"]
    manifest = manifest[[c for c in man_cols if c in manifest.columns]].sort_values(
        ["included", "track", "experiment_accession"], ascending=[False, True, True])
    manifest.to_csv(os.path.join(paths.out, "encode_vtrna1_source_manifest.csv"), index=False)

    # gene-ID audit
    audit_rows = []
    gtf_df = pd.DataFrame(gtf_records)
    for (track, asm, ann, schema), grp in wide.groupby(["track", "genome_assembly", "annotation_version", "quant_schema"], dropna=False):
        for sym in TARGET_SYMBOLS:
            p = sym.replace("-", "_")
            stat = grp[f"{p}_measurement_status"]
            present = stat.isin(["present", "present_effective_length_0"])
            reported = sorted(set(grp.loc[present, f"{p}_gene_id"].dropna()))
            gtf_ids = sorted(set(gtf_df[(gtf_df.annotation_label == f"{asm} {ann}") & (gtf_df.gene_symbol == sym)]["gene_id_gtf"])) if not gtf_df.empty else []
            audit_rows.append({
                "track": track, "pipeline_schema": schema,
                "analysis_titles": join(grp["analysis_title"].dropna().unique()),
                "annotation_version_source": join(grp["annotation_version_source"].dropna().unique()),
                "genome_assembly": asm, "annotation_version": ann,
                "reference_gtf": join(gtf_df[gtf_df.annotation_label == f"{asm} {ann}"]["gtf_accession"].unique()) if not gtf_df.empty else "",
                "gene_symbol": sym,
                "gene_id_in_annotation_gtf": join(gtf_ids),
                "gene_id_as_reported": join(reported),
                "stable_gene_id": join(sorted({strip_version(x) for x in reported + gtf_ids})),
                "number_of_files": len(grp),
                "number_of_files_present": int(present.sum()),
                "number_of_files_absent": int((stat.str.startswith("NA_absent")).sum()),
                "number_of_files_unavailable_or_unparsed": int(stat.isin(["NA_file_unavailable", "NA_parse_failed"]).sum()),
                "number_present_with_effective_length_0": int((stat == "present_effective_length_0").sum()),
            })
    # annotations encountered only among non-selected (alternate) files, e.g. hg19 V19 / GRCh38 V24
    alt = manifest[manifest.row_type == "alternate_file_not_selected"] if not manifest.empty else manifest
    seen_keys = {(r["track"], r["genome_assembly"], r["annotation_version"]) for r in audit_rows}
    if len(alt):
        for (track, asm, ann), grp in alt.groupby(["track", "genome_assembly", "annotation_version"], dropna=False):
            if (track, asm, ann) in seen_keys:
                continue
            lab = f"{asm} {ann}"
            for sym in TARGET_SYMBOLS:
                g = gtf_df[(gtf_df.annotation_label == lab) & (gtf_df.gene_symbol == sym)] if not gtf_df.empty else gtf_df
                audit_rows.append({
                    "track": track, "pipeline_schema": "not parsed (alternate files not selected)",
                    "analysis_titles": join(grp["analysis_title"].dropna().unique()),
                    "annotation_version_source": "file.genome_annotation", "genome_assembly": asm,
                    "annotation_version": ann,
                    "reference_gtf": join(g["gtf_accession"].unique()) if len(g) else "",
                    "gene_symbol": sym, "gene_id_in_annotation_gtf": join(g["gene_id_gtf"]) if len(g) else "",
                    "gene_id_as_reported": "", "stable_gene_id": join(g["stable_gene_id"].unique()) if len(g) else "",
                    "number_of_files": 0, "number_of_files_present": 0, "number_of_files_absent": 0,
                    "number_of_files_unavailable_or_unparsed": 0, "number_present_with_effective_length_0": 0,
                    "n_alternate_files_not_selected": int(len(grp)),
                })
    audit_df = pd.DataFrame(audit_rows)
    audit_df.to_csv(os.path.join(paths.out, "vtrna1_gene_identifier_audit.csv"), index=False)

    # QC + summaries
    qc_lines = qc_and_summaries(paths, md, experiments, classified, manifest, tot, sml, results, failed, audit_df)
    meta_info = {"run_started_utc": run_started, "run_finished_utc": utcnow(),
                 "encode_base": BASE, "assay_titles_considered": ASSAY_TITLES_CONSIDERED,
                 "n_experiments_queried": len(experiments)}
    with open(os.path.join(paths.meta, "run_info.json"), "w") as fh:
        json.dump(meta_info, fh, indent=2)
    print("\n".join(qc_lines))


# ---------------------------------------------------------------------------
# QC and summaries
# ---------------------------------------------------------------------------


def gene_summary(df: pd.DataFrame, value: str, track: str) -> pd.DataFrame:
    rows = []
    for sym in TARGET_SYMBOLS:
        p = sym.replace("-", "_")
        stat = df[f"{p}_measurement_status"]
        vals = pd.to_numeric(df[f"{p}_{value}"], errors="coerce")
        measured = vals.notna()
        quant = (stat == "present") & measured
        v = vals[measured]
        vq = vals[quant & measured]
        imax = v.idxmax() if len(v) and v.max() > 0 else None
        rows.append({
            "track": track, "quant_schema": join(df["quant_schema"].dropna().unique()), "gene_symbol": sym, "metric": value,
            "n_samples": len(df), "n_measured": int(measured.sum()),
            "n_quantifiable": int(quant.sum()),
            f"n_{value}_gt_0": int((v > 0).sum()),
            f"median_{value}": float(v.median()) if len(v) else np.nan,
            f"mean_{value}": float(v.mean()) if len(v) else np.nan,
            f"max_{value}": float(v.max()) if len(v) else np.nan,
            f"median_{value}_quantifiable_only": float(vq.median()) if len(vq) else np.nan,
            "biosample_with_max": df.loc[imax, "biosample_term_name"] if imax is not None else "",
            "file_with_max": df.loc[imax, "file_accession"] if imax is not None else "",
        })
    return pd.DataFrame(rows)


def biosample_summary(df: pd.DataFrame, value: str) -> pd.DataFrame:
    agg = {"N_biological_samples": ("file_accession", "count"),
           "N_experiments": ("experiment_accession", "nunique"),
           "N_donors": ("donor_accession", lambda s: len({x for v in s.dropna() for x in str(v).split(";") if x}))}
    for sym in TARGET_SYMBOLS:
        p = sym.replace("-", "_")
        agg[f"median_{value}_{p}"] = (f"{p}_{value}", lambda s: pd.to_numeric(s, errors="coerce").median())
        agg[f"n_{value}_gt0_{p}"] = (f"{p}_{value}", lambda s: int((pd.to_numeric(s, errors="coerce") > 0).sum()))
    out = df.groupby(["biosample_type", "biosample_term_name"], dropna=False).agg(**agg).reset_index()
    return out.sort_values(["biosample_type", "biosample_term_name"])


def qc_and_summaries(paths, md, experiments, classified, manifest, tot, sml, results, failed, audit_df) -> list[str]:
    L = []
    P = L.append
    P("=" * 78)
    P("QC REPORT - ENCODE primary human VTRNA1 atlas")
    P("=" * 78)
    P(f"Run date (UTC): {utcnow()}")
    tracks = Counter(classified[e['@id']]["track"] for e in experiments)
    P(f"Total ENCODE experiments queried (human, released, RNA-seq family): {len(experiments)}")
    P("  by assay_title: " + ", ".join(f"{k}={v}" for k, v in Counter(e.get('assay_title') for e in experiments).most_common()))
    P(f"Experiments passing all inclusion criteria (any track): {sum(c['eligible'] for c in classified.values())}")
    for track, label, df in [("total_RNA", "total RNA-seq (PRIMARY)", tot), ("small_RNA", "small RNA-seq (SEPARATE)", sml)]:
        P("")
        P(f"--- Track: {label} ---")
        elig = [e for e in experiments if classified[e['@id']]["track"] == track and classified[e['@id']]["eligible"]]
        P(f"Experiments with this assay: {tracks.get(track, 0)}")
        P(f"Experiments passing biosample/perturbation criteria: {len(elig)}")
        m = manifest[manifest.track == track]
        no_q = m[m.row_type == "eligible_no_gene_quantification"]["experiment_accession"].nunique()
        P(f"  of which with no usable released gene quantification: {no_q}")
        P(f"Experiments contributing quantifications: {df['experiment_accession'].nunique() if len(df) else 0}")
        P(f"Unique biosample terms: {df['biosample_term_name'].nunique() if len(df) else 0}")
        P(f"Unique donors: {len({x for v in df['donor_accession'].dropna() for x in str(v).split(';') if x}) if len(df) else 0}")
        P(f"Biological samples/replicates (rows): {len(df)}")
        P(f"Quantification files selected/downloaded OK/parsed OK: {len(df)}/"
          f"{int(df['download_ok'].fillna(False).sum()) if len(df) else 0}/{int(df['parsed_ok'].fillna(False).sum()) if len(df) else 0}")
        if not len(df):
            continue
        P("Assembly/annotation/schema: " + "; ".join(
            f"{a} {b} [{c}]: {n}" for (a, b, c), n in df.groupby(['genome_assembly', 'annotation_version', 'quant_schema'], dropna=False).size().items()))
        P("Biosample types: " + ", ".join(f"{k}={v}" for k, v in df['biosample_type'].value_counts().items()))
        for sym in TARGET_SYMBOLS:
            p = sym.replace("-", "_")
            st = df[f"{p}_measurement_status"]
            tpm = pd.to_numeric(df[f"{p}_TPM"], errors="coerce")
            cnt = pd.to_numeric(df[f"{p}_count_sense"], errors="coerce")
            P(f"  {sym}: row present in {int(st.isin(['present', 'present_effective_length_0']).sum())}/{len(df)} files; "
              f"status {st.value_counts().to_dict()}; TPM>0 in {int((tpm > 0).sum())} (of {int(tpm.notna().sum())} with TPM); "
              f"sense read count>0 in {int((cnt > 0).sum())} (of {int(cnt.notna().sum())} with counts)")
        missing = df[df[[f"{s.replace('-', '_')}_measurement_status" for s in TARGET_SYMBOLS]].apply(lambda r: r.astype(str).str.startswith("NA_").any(), axis=1)]
        P(f"Files with any VTRNA1 gene missing/NA: {len(missing)}" + (" -> " + ", ".join(missing['file_accession'].head(40)) if len(missing) else ""))
        P(f"Duplicate checks: duplicate file rows={df['file_accession'].duplicated().sum()}; "
          f"duplicate (experiment, bio-rep)={df.duplicated(['experiment_accession', 'biological_replicate']).sum()}; "
          f"library reused across rows={df['library_accession'].duplicated().sum()}; "
          f"biosample reused across rows={df['biosample_accession'].duplicated().sum()}")
        if df['biosample_accession'].duplicated().any():
            d = df[df['biosample_accession'].duplicated(keep=False)].groupby('biosample_accession')['experiment_accession'].apply(lambda s: ';'.join(sorted(set(s))))
            P(f"  biosamples appearing in >1 row (same material in separate experiments/replicates): {len(d)}")
            for k, v in list(d.items())[:20]:
                P(f"    {k}: {v}")
        g = df.groupby(["donor_accession", "biosample_term_name"]).agg(n=("file_accession", "size"), n_exp=("experiment_accession", "nunique"), n_lab=("lab", "nunique"))
        P(f"Donor x biosample-term groups with >1 row (distinct ENCODE biosamples from the same donor & cell type/tissue; "
          f"NOT independent donors): {(g.n > 1).sum()} groups covering {int(g.n[g.n > 1].sum())} rows; spanning >1 experiment: {(g.n_exp > 1).sum()}; spanning >1 lab: {(g.n_lab > 1).sum()}; "
          f"largest: " + ", ".join(f"{a}/{b}={n}" for (a, b), n in g.n.sort_values(ascending=False).head(4).items()))
        P(f"Rows whose selected file pools >1 technical replicate: {df['file_technical_replicates'].astype(str).str.contains(';').sum()}")
        P(f"Rows where other technical-replicate files of the same biological replicate exist but were NOT used "
          f"(one file kept per biological replicate): {df['other_technical_replicate_files_not_selected'].fillna('').astype(bool).sum()}")
        rs = df[df.quant_schema == "RSEM_genes_results"]
        if len(rs):
            s_ = pd.to_numeric(rs['file_sum_TPM'], errors='coerce')
            P(f"[RSEM files n={len(rs)}] per-file sum(TPM) min={s_.min():.0f} median={s_.median():.0f} max={s_.max():.0f} (expect ~1e6); "
              f"genes/file {int(rs['file_n_genes'].min())}-{int(rs['file_n_genes'].max())}")
            for sym in TARGET_SYMBOLS:
                p = sym.replace("-", "_")
                el = pd.to_numeric(rs[f"{p}_effective_length"], errors="coerce")
                ln = pd.to_numeric(rs[f"{p}_length"], errors="coerce")
                t = pd.to_numeric(rs[f"{p}_TPM"], errors="coerce")
                c = pd.to_numeric(rs[f"{p}_expected_count"], errors="coerce")
                P(f"  {sym}: length={sorted(ln.dropna().unique().tolist())}; effective_length>0 in {int((el > 0).sum())}/{int(el.notna().sum())}; "
                  f"TPM q(0,50,90,99,100)={np.round(np.nanpercentile(t, [0, 50, 90, 99, 100]), 3).tolist() if t.notna().any() else 'NA'}; "
                  f"expected_count q={np.round(np.nanpercentile(c, [0, 50, 90, 99, 100]), 2).tolist() if c.notna().any() else 'NA'}; "
                  f"negative={int(((t < 0) | (c < 0)).sum())}")
        st_ = df[df.quant_schema == "STAR_ReadsPerGene"]
        if len(st_):
            s_ = pd.to_numeric(st_['file_total_assigned_counts'], errors='coerce')
            P(f"[STAR count files n={len(st_)}] assigned gene counts/file min={s_.min():.0f} median={s_.median():.0f} max={s_.max():.0f}; "
              "sense column: " + ", ".join(f"{k}={v}" for k, v in st_['file_sense_count_column'].value_counts().items()))
            for sym in TARGET_SYMBOLS:
                p = sym.replace("-", "_")
                c = pd.to_numeric(st_[f"{p}_count_sense"], errors="coerce")
                cp = pd.to_numeric(st_[f"{p}_CPM"], errors="coerce")
                P(f"  {sym}: sense count q(0,50,90,99,100)={np.round(np.nanpercentile(c, [0, 50, 90, 99, 100]), 1).tolist() if c.notna().any() else 'NA'}; "
                  f"CPM q={np.round(np.nanpercentile(cp, [0, 50, 90, 99, 100]), 3).tolist() if cp.notna().any() else 'NA'}; negative={int((c < 0).sum())}")
        P("Library classes: " + "; ".join(f"{k}={v}" for k, v in df['library_rna_class'].value_counts().items()))
        P("Vault-RNA suitability: " + "; ".join(f"{k}={v}" for k, v in df['vault_rna_suitability'].value_counts().items()))
    P("")
    P(f"Failed downloads: {len(failed)}" + ("" if not failed else " -> " + ", ".join(a for a, _ in failed)))
    P("")
    P("--- Gene identifier audit ---")
    with pd.option_context("display.width", 300, "display.max_columns", 50, "display.max_colwidth", 40):
        P(audit_df.to_string(index=False))

    def emit_block(df, title, value, gene_file, bs_file):
        P("")
        P("=" * 78)
        P(title)
        P("=" * 78)
        if not len(df):
            P("(no rows)")
            return
        gs = gene_summary(df, value, df["track"].iloc[0])
        gs.to_csv(os.path.join(paths.out, gene_file), index=False)
        with pd.option_context("display.width", 300, "display.max_columns", 50):
            P(gs.drop(columns=["track"]).to_string(index=False))
        bsum = biosample_summary(df, value)
        bsum.to_csv(os.path.join(paths.out, bs_file), index=False)
        P("")
        P(f"Tissue / cell-type summary (median {value}):")
        with pd.option_context("display.max_rows", 1000, "display.width", 300):
            P(bsum[["biosample_type", "biosample_term_name", "N_biological_samples"]
                   + [f"median_{value}_{s.replace('-', '_')}" for s in TARGET_SYMBOLS]].to_string(index=False))

    emit_block(tot[tot.quant_schema == "RSEM_genes_results"],
               "SUMMARY A: total RNA-seq, ENCODE uniform pipeline (STAR+RSEM) - TPM",
               "TPM", "summary_vtrna1_by_gene_totalRNA_RSEM_TPM.csv", "summary_vtrna1_by_biosample_totalRNA_RSEM_TPM.csv")
    emit_block(tot[tot.quant_schema == "STAR_ReadsPerGene"],
               "SUMMARY B: total RNA-seq, UW/Altius lab pipeline (STAR raw counts, GENCODE V43) - CPM (derived)",
               "CPM", "summary_vtrna1_by_gene_totalRNA_STARcounts_CPM.csv", "summary_vtrna1_by_biosample_totalRNA_STARcounts_CPM.csv")
    emit_block(sml, "SUMMARY C: small RNA-seq (<200 nt), ENCODE small-RNA pipeline (STAR counts) - CPM (derived)",
               "CPM", "summary_vtrna1_by_gene_smallRNA_CPM.csv", "summary_vtrna1_by_biosample_smallRNA_CPM.csv")
    with open(os.path.join(paths.out, "qc_report.txt"), "w") as fh:
        fh.write("\n".join(L) + "\n")
    return L


if __name__ == "__main__":
    main()
