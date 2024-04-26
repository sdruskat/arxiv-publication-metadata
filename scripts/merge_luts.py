import json
import logging

merged_lut = dict()

log = logging.getLogger(__name__)
file_handler = logging.FileHandler(f"{snakemake.output[0]}.log", mode="w")
file_handler.setFormatter(logging.Formatter(
    fmt="[%(asctime)s] [%(levelname)8s] --- %(message)s (%(module)s.%(funcName)s > %(filename)s:%(lineno)s)"
))
log.setLevel(logging.getLevelName("DEBUG"))
log.addHandler(file_handler)

total = 0
_sum = 0

def deep_merge_lut(exist: dict, new: dict) -> dict:
    global _sum, total
    _sum += len(new)
    exist = exist | new
    total = len(exist)
    return exist | new


for lut_file in snakemake.input:
    with open(lut_file, "r") as fi:
        data = json.load(fi)

    merged_lut = deep_merge_lut(merged_lut, data)

print(f"Total number of items = {total} vs. sum of single LUTs = {_sum}")
size = merged_lut.__sizeof__()
print(f"LUT has size {size} bytes ({size / 1e+6} MB, {size / 1e+9} GB).")

with open(snakemake.output[0], "w") as of:
    json.dump(merged_lut, of)
