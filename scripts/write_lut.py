# SPDX-FileCopyrightText: 2024 German Aerospace Center (DLR)
# SPDX-FileContributor: Stephan Druskat <stephan.druskat@dlr.de>
#
# SPDX-License-Identifier: MIT

import json
import xml.etree.ElementTree as ET
from datetime import datetime
import logging

import xmlschema


log = logging.getLogger(__name__)

file_handler = logging.FileHandler(f"{snakemake.output[0]}.log", mode="w")
file_handler.setFormatter(logging.Formatter(
    fmt="[%(asctime)s] [%(levelname)8s] --- %(message)s (%(module)s.%(funcName)s > %(filename)s:%(lineno)s)"
))
log.setLevel(logging.getLevelName("DEBUG"))
log.addHandler(file_handler)

def get_lut_for_file():
    xml_file = snakemake.input.xml
    arxiv_xsd_file = snakemake.input.arxiv_xsd
    xsd = xmlschema.XMLSchema(arxiv_xsd_file)

    prefix = "{http://arxiv.org/OAI/arXivRaw/}"

    tree = ET.parse(xml_file)

    _lut = dict()

    for record in tree.find("ListRecords"):
        metadata = record.find("{http://www.openarchives.org/OAI/2.0/}metadata")
        if metadata:
            arxiv_data = metadata.find(f"{prefix}arXivRaw")
            data, errors = xsd.to_dict(arxiv_data, validation="lax")
            if len(data) > 0:
                arxiv_id = data[f"{prefix}id"]
                versions = data[f"{prefix}version"]
                if not arxiv_id or not versions:
                    log.warning(f"Metadata did not contain both versions and ID in {xml_file}.")
                else:
                    for version in versions:
                        v = version["@version"]
                        date = version[f"{prefix}date"]
                        if date:
                            date_obj = datetime.strptime(date, "%a, %d %b %Y %H:%M:%S %Z")
                            date_f = datetime.strftime(date_obj, "%Y-%m-%d")
                            _lut[f"{arxiv_id}{v}"] = date_f
                            log.debug(f"Recorded {date_f} for {arxiv_id}{v}")

    return _lut

if __name__ == '__main__':
    lut_for_file = get_lut_for_file()
    with open(snakemake.output.lut, "w") as out:
            json.dump(lut_for_file, out)
    if not lut_for_file:
        log.warning(f"LUT for {snakemake.input.xml} was empty.")
