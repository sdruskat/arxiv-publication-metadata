import logging
import json
from datetime import datetime
import xml.etree.ElementTree as ET

import requests
import xmlschema
from ratelimit import limits, sleep_and_retry


log = logging.getLogger(__name__)

file_handler = logging.FileHandler(f"{snakemake.output[0]}.log", mode="w")
file_handler.setFormatter(logging.Formatter(
    fmt="[%(asctime)s] [%(levelname)8s] --- %(message)s (%(module)s.%(funcName)s > %(filename)s:%(lineno)s)"
))
log.setLevel(logging.getLevelName("DEBUG"))
log.addHandler(file_handler)


def extract_lut_from_xml(xml_source: str, xsd_file: str, is_file: bool = True):
    xsd = xmlschema.XMLSchema(xsd_file)

    prefix = "{http://arxiv.org/OAI/arXivRaw/}"
    oai_prefix = "{http://www.openarchives.org/OAI/2.0/}"


    if is_file:
        tree = ET.parse(xml_source)
    else:
        tree = ET.fromstring(xml_source)

    _lut = dict()

    record = tree.find(f"{oai_prefix}GetRecord").find(f"{oai_prefix}record")
    metadata = record.find(f"{oai_prefix}metadata")

    if metadata:
        arxiv_data = metadata.find(f"{prefix}arXivRaw")
        data, errors = xsd.to_dict(arxiv_data, validation="lax")
        if data:
            arxiv_id = data[f"{prefix}id"]
            versions = data[f"{prefix}version"]
            if not arxiv_id or not versions:
                log.warning(f"Metadata did not contain both versions and ID in {xml_source}.")
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

@sleep_and_retry
@limits(calls=1, period=180)
def assert_versions() -> list[str]:
    """
    Test if all IDs in all Extract-URLs ArXiv JSONs are in the compiled LUT.
    Return missing IDs.

    Must only run once every 3 seconds.
    """

    with open(snakemake.input.lut, "r") as lut:
        lut_data = json.load(lut)

    missing = []

    for jf in snakemake.input.arxiv_urls:
        file_ym = jf.split("/")[-1].split(".")[0]
        with open(jf, "r") as ji:
            url_data = json.load(ji)
        ym_data = url_data[file_ym]
        files = ym_data["files"]
        for pdf_name in files.keys():
            ver_id = pdf_name.split(".pdf")[0]
            if not ver_id in lut_data.keys():
                missing.append(ver_id)

    for missing_ver_id in missing:
        split_id = missing_ver_id.split("v")
        record_id = split_id[0]
        ver_id = split_id[-1]
        _ver_id_int = int(ver_id)
        response = requests.get(f"http://export.arxiv.org/oai2?verb=GetRecord&identifier=oai:arXiv.org:{record_id}&metadataPrefix=arXivRaw")
        if response.status_code == 200:
            lut_from_response = extract_lut_from_xml(response.text, snakemake.input.arxiv_xsd, False)
            lut_data = lut_data | lut_from_response

    return lut_data

if __name__ == '__main__':
    patched_lut = assert_versions()
    with open(snakemake.output[0], "w") as out:
        json.dump(patched_lut, out)