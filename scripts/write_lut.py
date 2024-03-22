import json
import xml.etree.ElementTree as ET
from xml.etree.ElementTree import Element as E
from datetime import datetime

NAMESPACES = {"oai": "http://www.openarchives.org/OAI/2.0/",
              "arxiv": "http://arxiv.org/OAI/arXivRaw/"}

def find(node: E, search_str: str, namespace: str = "oai"):
    return node.find(f"{namespace}:{search_str}", NAMESPACES)

def findall(node: E, search_str: str, namespace: str = "oai"):
    return node.findall(f"{namespace}:{search_str}", NAMESPACES)


tree = ET.parse(snakemake.input[0])
root = tree.getroot()

lut = dict()

list_records = root.find("ListRecords")
for record in findall(list_records, "record"):
    metadata = find(record, "metadata")
    arxiv = find(metadata, "arXivRaw", "arxiv")
    arxiv_id = find(arxiv, "id", "arxiv")
    versions = findall(arxiv, "version", "arxiv")
    for version in versions:
        v = version.attrib["version"]
        date = find(version, "date", "arxiv")
        date_obj = datetime.strptime(date.text, "%a, %d %b %Y %H:%M:%S %Z")
        date_f = datetime.strftime(date_obj, "%Y-%m-%d")
        lut[f"{arxiv_id.text}{v}"] = date_f

with open(snakemake.output[0], "w") as of:
    json.dump(lut, of)