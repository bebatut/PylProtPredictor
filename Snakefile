from snakemake.remote.FTP import RemoteProvider as FTPRemoteProvider
from snakemake.utils import available_cpu_count

FTP = FTPRemoteProvider()

configfile: "config.yaml"


rule help:
    '''Print the help message'''
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
    '''Print the help message'''
    shell: "echo 0.2"


rule purge:
    '''Delete the output directory'''
    input:
        expand(
            "{output_dir}",
            output_dir=config["output_dir"])
    message: "!!! DELETING THE OUTPUT DIRECTORY '{input}' !!!"
    shell: "rm -rf {input}"


rule all:
    '''Run the complete pipeline'''
    input:
        expand(
            "{output_dir}/{file}",
            output_dir=config["output_dir"],
            file="report.html")


rule PylProtPredictor:
    '''Run the complete pipeline'''
    input: rules.all.input


rule prepare_database:
    '''Download and format the Uniref90 database for Diamond'''
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
    '''Predict CDS from the input genome using Prodigal'''
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
    '''Predict potential PYL-contrainin proteins from predicted CDS'''
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
        percentage_info=expand(
            "{output_dir}/{file}",
            output_dir=config["output_dir"],
            file="percentage_info.csv"),
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
    '''Align predicted Pyl proteins on Uniref90 using Diamond'''
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
    '''Validate potential Pyl protein by analysing the alignments'''
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
    '''Print a quick HTML report'''
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
