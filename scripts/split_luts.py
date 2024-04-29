import json
import logging
from pathlib import Path

log = logging.getLogger(__name__)

file_handler = logging.FileHandler(f"{snakemake.output[0]}.log", mode="w")
file_handler.setFormatter(logging.Formatter(
    fmt="[%(asctime)s] [%(levelname)8s] --- %(message)s (%(module)s.%(funcName)s > %(filename)s:%(lineno)s)"
))
log.setLevel(logging.getLevelName("DEBUG"))
log.addHandler(file_handler)


patched_lut = snakemake.input[0]
out_dir = snakemake.output[0]

luts = dict()

with open(patched_lut, "r") as lutf:
    data = json.load(lutf)

    data_len = len(data)

    for i, (k, v) in enumerate(data.items()):
        if i % 100000 == 0:
            log.debug(f"Splitting item {i} of {data_len} ({i / data_len * 100}%)")
        file_name = k.split("/")[0].split(".")[0] + ".json"
        if file_name in luts:
            luts[file_name] = luts[file_name] | {k: v}
        else:
            luts[file_name] = {k: v}

luts_count = len(luts)

od = Path(out_dir)
od.mkdir(parents=True, exist_ok=True)

for i, (file, lut) in enumerate(luts.items()):
    log.debug(f"Writing LUT {i + 1} ({file}) of {luts_count}")
    with open(f"{out_dir}/{file}", "w") as out:
        json.dump(lut, out)