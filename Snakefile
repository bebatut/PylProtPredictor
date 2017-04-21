from snakemake.remote.FTP import RemoteProvider as FTPRemoteProvider
FTP = FTPRemoteProvider()

configfile: "config.yaml"

rule all:
    input:
        expand("{dirpath}/{cds}", dirpath=config["dirpath"], cds="predicted_cds.fasta")

rule download_uniref90:
    input:
        # only keeping the file so we can move it out to the cwd
        FTP.remote("ftp.uniprot.org/pub/databases/uniprot/uniref/uniref90/uniref90.fasta.gz", keep_local=True)
    output:
        "data/uniref90.fasta"
    shell:
        "gunzip {input} | mv {output}"

rule predict_cds:
    input:
        expand("{dirpath}/{genome}", dirpath=config["dirpath"], genome=config["genome"])
    output:
        cds=expand("{dirpath}/{cds}", dirpath=config["dirpath"], cds="predicted_cds.fasta"),
        info=expand("{dirpath}/{info}", dirpath=config["dirpath"], info="cds_prediction_info.txt"),
        log=expand("{dirpath}/{log}", dirpath=config["dirpath"], log="cds_prediction_log.txt")
    shell:
        "prodigal"
        " -i {input}"
        " -d {output.cds}"
        " -f gbk"
        " -g 11"
        " -o {output.info}"
        " 2> {output.log}"
