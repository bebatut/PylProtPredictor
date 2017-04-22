from snakemake.remote.FTP import RemoteProvider as FTPRemoteProvider
FTP = FTPRemoteProvider()

configfile: "config.yaml"


rule all:
    input:
        expand(
            "{output_dir}/{cds}",
            output_dir=config["output_dir"],
            cds="potential_pyl_sequences.fasta")


rule download_uniref90:
    input:
        # only keeping the file so we can move it out to the cwd
        FTP.remote(
            config["ref_database_url"],
            keep_local=True)
    output:
        expand(
            "{data_dir}/{ref_database}.fasta",
            data_dir=config["data_dir"],
            ref_database=config["ref_database"])
    shell:
        "gunzip {input} | mv {output}"


rule prepare_uniref90:
    input:
        expand(
            "{data_dir}/{ref_database}.fasta",
            data_dir=config["data_dir"],
            ref_database=config["ref_database"])
    output:
        expand(
            "{data_dir}/{ref_database}.dmnd",
            data_dir=config["data_dir"],
            ref_database=config["ref_database"])
    shell:
        "diamond makedb"
        " --in {input}"
        " --db {output}"
        " --quiet"


rule predict_cds:
    input:
        config["genome"]
    output:
        cds=expand(
            "{output_dir}/{cds}",
            output_dir=config["output_dir"],
            cds="predicted_cds.fasta"),
        info=expand(
            "{output_dir}/{info}",
            output_dir=config["output_dir"],
            info="cds_prediction_info.txt"),
        log=expand(
            "{output_dir}/{log}",
            output_dir=config["output_dir"],
            log="cds_prediction_log.txt")
    shell:
        "prodigal"
        " -i {input}"
        " -d {output.cds}"
        " -f gbk"
        " -g 11"
        " -o {output.info}"
        " 2> {output.log}"


rule predict_potential_pyl_proteins:
    input:
        genome=config["genome"],
        predicted_cds=expand(
            "{output_dir}/{cds}",
            output_dir=config["output_dir"],
            cds="predicted_cds.fasta")
    output:
        potential_pyl_seq=expand(
            "{output_dir}/{pot_pyl_seq}",
            output_dir=config["output_dir"],
            pot_pyl_seq="potential_pyl_sequences.fasta"),
        log=expand(
            "{output_dir}/{log}",
            output_dir=config["output_dir"],
            log="pyl_protein_prediction_log.txt")
    script:
        "src/predict_pyl_proteins.py"


rule report:
    input:
        ""
    output:
        "report.html"
    run:
        from snakemake.utils import report
        with open(input[0]) as vcf:
            n_calls = sum(1 for l in vcf if not l.startswith("#"))

        report("""
        Prediction of PYL proteins
        ==========================

        This resulted in {n_calls} variants (see Table T1_).
        """, output[0], T1=input[0])
