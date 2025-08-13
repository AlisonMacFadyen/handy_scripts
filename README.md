# Handy Scripts

A collection of small but useful scripts I’ve written to solve research and data-related problems.

They cover **comparative genomics, data manipulation, file processing for uploading to data repositories and protein-related bioinformatics tasks**.

This repository serves as both a personal toolbox and a showcase of my coding style, documentation and problem-solving approach.

---

## Repository Structure

| Directory                  | Scripts                                                                          | Description                                                                                                                                                                                    |
| -------------------------- | -------------------------------------------------------------------------------- | ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| **comparative\_genomics/** | `process_ariba_phandango.py`, `process_snippy.py`                                | Scripts for preparing and transforming comparative genomics output from tools like **ARIBA** and **Snippy** for downstream analysis or visualisation.                          |
| **data\_manipulation/**    | `output_unique.py`                                                               | Utility to extract unique values from datasets, useful for deduplication, quick data filtering and list processing.                                                                          |
| **data\_processing/**      | `prepare_ena_manifest.py`, `process_read_paths.py`, `upload_to_ENA.sh`           | Tools for managing sequence data submission workflows, including generating **ENA** manifest files, processing read file paths and automating uploads to the **European Nucleotide Archive**. |
| **protein\_associated/**   | `combine_aa_AF2_one_vs_one.py`, `extract_domains.py`, `parse_hmmscan_results.py` | Protein-focused bioinformatics scripts for combining amino acid comparison outputs (AlphaFold2 one-vs-one runs), extracting protein domains and parsing HMMER `hmmscan` output.               |

---

## Usage

Each script is self-contained. You can run them directly with Python:

```bash
python path/to/script.py
```

For shell scripts:

```bash
bash path/to/script.sh
```

Some scripts may require dependencies — check the script’s import lines for details.

---

## Skills Demonstrated

* **Bioinformatics workflows** (genomics, protein analysis, ENA submissions)
* **Data parsing and transformation** (CSV, TSV, FASTA and tool-specific formats)
* **Automation** of repetitive research tasks
* **Filesystem operations** and batch processing
* **Use of external bioinformatics tools** such as ARIBA, Snippy, HMMER and AlphaFold2

---

## License

This repository is licensed under the terms of the MIT license, feel free to use and adapt these scripts.
