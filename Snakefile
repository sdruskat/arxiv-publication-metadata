rule all:
    input:
        expand("data/extracted/{file_name}.xml", file_name=glob_wildcards(os.path.join("data", "{file_name}.xml.gz")).file_name)

rule gunzip:
    input:
        "data/{file_name}.xml.gz"
    output:
        "data/extracted/{file_name}.xml"
    run:
        shell("mkdir -p data/extracted && gunzip -c {input} > {output}")

rule clean:
    threads: 1
    shell:
        "rm -rf data/extracted .cache/ .conda/ .snakemake/ logs/"