from snakemake.remote.FTP import RemoteProvider as FTPRemoteProvider
from snakemake.utils import available_cpu_count

FTP = FTPRemoteProvider()

configfile: "config.yaml"

rule help:
    threads: available_cpu_count()
    message:
      "------------------------------------------------" +"\n"
      "    Welcome to the PylProtPredictor pipeline    " +"\n"
      "------------------------------------------------" +"\n"
      # "" +"\n"
      "  version:            0.2" +"\n"
      "  max_threads:        {threads}" +"\n"
      "  database location:  "+config["data_dir"]+"/"+config["ref_database"]+".dmnd" +"\n"
      "  input genome        "+config["genome"] +"\n"
      "  output directory:   "+config["output_dir"]+"/" +"\n"
      "" +"\n"
      "Usage:" +"\n"
      " - Print this help message:" +"\n"
      "      snakemake help" +"\n"
      " - List available rules:" +"\n"
      "      snakemake --list" +"\n"
      " - Run the complete pipeline:" +"\n"
      "      snakemake --cores -p all" +"\n"
      " - Only prepare the database:" +"\n"
      "      snakemake --cores -p prepare_database" +"\n"
      " - Predict potential Pyl proteins (does not check predictions):" +"\n"
      "      snakemake --cores -p predict_potential_pyl_proteins" +"\n"
      "------------------------------------------------" +"\n"


rule version:
    shell: "echo 0.2"


rule purge:
    input:
        expand(
            "{output_dir}",
            output_dir=config["output_dir"])
    message: "!!! DELETING THE OUTPUT DIRECTORY '{input}' !!!"
    shell: "rm -rf {input}"


rule all:
    input:
        expand(
            "{output_dir}/{file}",
            output_dir=config["output_dir"],
            file="report.html")


rule PylProtPredictor:
    """alias for 'all' rule"""
    input: rules.all.input


rule prepare_database:
    input:
        FTP.remote(
            config["ref_database_url"],
            keep_local=False)
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
            "{output_dir}/{file}",
            output_dir=config["output_dir"],
            file="predicted_cds.fasta"),
        info=expand(
            "{output_dir}/{file}",
            output_dir=config["output_dir"],
            file="cds_prediction_info.txt"),
        log=expand(
            "{output_dir}/{file}",
            output_dir=config["output_dir"],
            file="cds_prediction_log.txt")
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
            "{output_dir}/{file}",
            output_dir=config["output_dir"],
            file="predicted_cds.fasta")
    output:
        potential_pyl_sequences=expand(
            "{output_dir}/{file}",
            output_dir=config["output_dir"],
            file="potential_pyl_sequences.fasta"),
        pyl_protein_prediction_log=expand(
            "{output_dir}/{file}",
            output_dir=config["output_dir"],
            file="pyl_protein_prediction_log.txt"),
        predicted_cds_info=expand(
            "{output_dir}/{file}",
            output_dir=config["output_dir"],
            file="predicted_cds.csv"),
        tag_ending_cds_info=expand(
            "{output_dir}/{file}",
            output_dir=config["output_dir"],
            file="tag_ending_cds.csv"),
        potential_pyl_protein_info=expand(
            "{output_dir}/{file}",
            output_dir=config["output_dir"],
            file="potential_pyl_proteins.csv")
    script:
        "src/predict_pyl_proteins.py"


rule search_similarity:
    threads: available_cpu_count()
    input:
        potential_pyl_seq=expand(
            "{output_dir}/{file}",
            output_dir=config["output_dir"],
            file="potential_pyl_sequences.fasta"),
        ref_db=expand(
            "{data_dir}/{ref_database}.dmnd",
            data_dir=config["data_dir"],
            ref_database=config["ref_database"])
    output:
        potential_pyl_similarity_search=expand(
            "{output_dir}/{file}",
            output_dir=config["output_dir"],
            file="potential_pyl_similarity_search.txt")
    shell:
        "diamond blastp"
        " -d {input.ref_db}"
        " -q {input.potential_pyl_seq}"
        " -o {output}"
        " -k 1"
        " -e 0.01"
        " -f 6"
        " -b 0.5"
        " -p {threads} "
        " --quiet"


rule check_pyl_proteins:
    input:
        potential_pyl_similarity_search=expand(
            "{output_dir}/{file}",
            output_dir=config["output_dir"],
            file="potential_pyl_similarity_search.txt"),
        potential_pyl_sequences=expand(
            "{output_dir}/{file}",
            output_dir=config["output_dir"],
            file="potential_pyl_sequences.fasta")
    output:
        conserved_potential_pyl_seq=expand(
            "{output_dir}/{file}",
            output_dir=config["output_dir"],
            file="conserved_potential_pyl_sequences.fasta"),
        conserved_potential_pyl_seq_info=expand(
            "{output_dir}/{file}",
            output_dir=config["output_dir"],
            file="conserved_potential_pyl_sequences.csv"),
        rejected_potential_pyl_seq_info=expand(
            "{output_dir}/{file}",
            output_dir=config["output_dir"],
            file="rejected_potential_pyl_sequences.csv")
    script:
        "src/check_pyl_proteins.py"


rule report:
    input:
        predicted_cds_info=expand(
            "{output_dir}/{file}",
            output_dir=config["output_dir"],
            file="predicted_cds.csv"),
        tag_ending_cds=expand(
            "{output_dir}/{file}",
            output_dir=config["output_dir"],
            file="tag_ending_cds.csv"),
        potential_pyl_proteins=expand(
            "{output_dir}/{file}",
            output_dir=config["output_dir"],
            file="potential_pyl_proteins.csv"),
        conserved_potential_pyl_sequences=expand(
            "{output_dir}/{file}",
            output_dir=config["output_dir"],
            file="conserved_potential_pyl_sequences.csv")
    output:
        report=expand(
            "{output_dir}/{file}",
            output_dir=config["output_dir"],
            file="report.html")
    script:
        "src/write_report.py"
