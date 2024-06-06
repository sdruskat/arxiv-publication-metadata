<!--
SPDX-FileCopyrightText: 2024 German Aerospace Center (DLR)
SPDX-FileContributor: Stephan Druskat <stephan.druskat@dlr.de>

SPDX-License-Identifier: CC0-1.0
-->

# Documentation

This workflow takes as input ArXiv publication OAI-PMH2 metadata in `arXivRaw` format harvested via 
[`metha`](https://github.com/miku/metha)[^1].

## Steps ("Rules")

1. The dataset [^2] is downloaded from Zenodo, and extracted to give access to the XML files containing 
the harvested metadata.
2. From each XML file, the following metadata is extracted and written into a lookup table that maps
the ArXiv version identifier (e.g., `2404.12345v1`) to the value for the version:
    - Publication date in `YYYY-MM-DD` format.
Each lookup table is saved to a JSON file.
3. The lookup tables are merged into a single table.
4. An attempt is made to patch missing metadata into the lookup table.
    - For all ArXiv version identifiers that are present in the *Extract-URLs* dataset [^3], it is checked if they are 
      in the lookup table, attempts are made to retrieve missing metadata from the single OAI-PMH record of
      the publication.
5. The patched lookup table produced in rule 4. is split into several files named `<identifier prefix>.json` where
`identifier prefix` is the `YYMM`-formatted prefix of the ArXiv publication identifier of each identifier with that 
prefix.
6. The files containing the LUT for a specific ArXiv identifier prefix are archived in a tarball that is saved in
`results/`, and all file names (without paths) as saved as a JSON array in `file_names.json`.

A graphical overview  of the rules is given below:

![Rulegraph of the rules described above, generated via `snakemake --rulegraph | dot -Tsvg > rulegraph.svg` run in the repository root.](../rulegraph.svg)

[^1]: Martin Czygan, Thomas Gersch, ACz-UniBi, Justin Kelly, Gunnar Þór Magnússon, dvglc, & Natanael Arndt. (2024). _miku/metha: metha 0.3.5_ (v0.3.5). Zenodo. doi:[10.5281/zenodo.11066532](https://doi.org/10.5281/zenodo.11066532)
[^2]: Druskat, S. (2024). _ArXiv OAI-PMH arXivRaw publication metadata_ (Version 1) [Data set]. Zenodo. doi:[10.5281/zenodo.11065282](https://doi.org/10.5281/zenodo.11065282)
[^3]: Escamilla, E. _Extract-URLs_ [Data set]. <https://github.com/elescamilla/Extract-URLs>