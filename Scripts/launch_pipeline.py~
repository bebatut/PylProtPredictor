#!/usr/bin/python2.7 
# -*- coding: utf-8 -*- 

import sys 
import os 
import re 
from Bio import SeqIO

##########
#METHODES#
##########

#Verifie si le fichier est au format fasta 
def verif_fasta(file): 
	fasta=False
	for record in SeqIO.parse(file,"fasta"):
		fasta=True		  
	return fasta

#Verifie si fichier contient plusieurs scaffolds 
def verif_scaffold(file): 
	input_file=open(file,"r") 
	count_seq=0
	for line in input_file: 
		if (re.search(r"^>",line)): 
			count_seq+=1
	if count_seq==1: 
		scaffold=False 
	else: 
		scaffold=True 
	return scaffold 	 

#GESTION DES FICHIERS#
genome_complet=sys.argv[1]
dir_result_prediction=sys.argv[2]
dir_result_blast=sys.argv[3]
db=sys.argv[4]
evalue=sys.argv[5] 
align=sys.argv[6]
cds_predits=dir_result_prediction+"/predicted_CDS.fasta"
CDS_TAG_plus=dir_result_prediction+"/CDS_TAG_normal.fasta" #Contient tous les CDS se terminant par TAG du brin+  
CDS_TAG_moins=dir_result_prediction+"/CDS_TAG_complement.fasta" #Contient tous les CDS se terminant par TAG du brin-
position_CDS_TAG_plus=dir_result_prediction+"/position_CDS_TAG_normal.txt" #Contient la position sur le genome des CDS TAG du brin+
position_CDS_TAG_moins=dir_result_prediction+"/position_CDS_TAG_complement.txt" #Contient la position sur le genome des CDS TAG du brin- 
fichier_longueur_scaffolds=dir_result_prediction+"/length_scaffolds.txt" 
position_cds_intergen_plus=dir_result_prediction+"/position_keep_regions_normal.txt" #Contient les positions des regions a conserver (CDS+region intergenique en aval) (brin+) 
position_cds_intergen_moins=dir_result_prediction+"/position_keep_regions_complement.txt" #Contient les positions des regions a conserver (CDS+region intergenique en amont) (brin-) 
region_cds_intergen_plus=dir_result_prediction+"/keep_regions_normal.fasta" #Contient les regions a conserver (brin+)  
region_cds_intergen_moins=dir_result_prediction+"/keep_regions_complement.fasta" #Contient les regions a conserver (brin-) 
region_cds_intergen_moins_RC=dir_result_prediction+"/keep_regions_complement_RC.fasta" #Regions reverse complementees pour le brin -
genes_possibles_plus=dir_result_prediction+"/potential_genes_normal.fasta" #Liste des genes possibles (brin+) 
genes_possibles_moins=dir_result_prediction+"/potential_genes_complement.fasta" #Liste des genes possibles (brin-) 
protein_no_pyl_plus=dir_result_prediction+"/potential_proteins_normal.fasta" #Genes possibles traduits (brin+) 
protein_no_pyl_moins=dir_result_prediction+"/potential_proteins_complement.fasta"#Genes possibles traduits (brin-) 
protein_pyl_plus=dir_result_prediction+"/potential_PYL_proteins_normal.fasta" #Genes possibles traduits avec les tag codants = O (brin+) 
protein_pyl_moins=dir_result_prediction+"/potential_PYL_proteins_complement.fasta" #Genes possibles traduits avec les tag codants = O (brin-)
blast_plus=dir_result_blast+"/blast_normal.out" #resultat brut blast +
blast_moins=dir_result_blast+"/blast_complement.out" #resultat brut blast - 
blast_plus_xml=dir_result_blast+"/blast_normal.xml" #resultat brut xml blast + (pr traitement BioPython) 
blast_moins_xml=dir_result_blast+"/blast_complement.xml" #resultat brut xml blast - 
traitement_blast_plus=dir_result_blast+"/sort_results_blast_normal.txt" 
traitement_blast_moins=dir_result_blast+"/sort_results_blast_complement.txt" 


#Gestion scripts 
home_directory=os.environ['HOME']
path_prodigal=home_directory+"/pyl_protein_prediction/Prodigal" #chemin vers dossier où se situe prodigal.linux
prediction="python "+home_directory+"/pyl_protein_prediction/Scripts/prediction_ORF_pyl_prodigal.py "+genome_complet+" "+cds_predits+" "+dir_result_prediction #lancement prediction genes pyl pr genome completement assemblé 
prediction_scaffold="python "+home_directory+"/pyl_protein_prediction/Scripts/prediction_ORF_pyl_prodigal_scaffold.py "+genome_complet+" "+cds_predits+" "+dir_result_prediction #lancement prediction genes pyl pr genome avc scaffolds 
blast="python "+home_directory+"/pyl_protein_prediction/Scripts/blast.py "+protein_pyl_plus+" "+protein_pyl_moins+" "+db+" "+dir_result_blast+" "+evalue+" "+align
#lancement blast (resultats bruts+simplifiés) 

#Verification fichiers 
fasta_genome=verif_fasta(genome_complet) 
scaffolds=verif_scaffold(genome_complet)
if not fasta_genome: 
	print "Complete genome : need fasta format" 
 
if fasta_genome:

	#1.Prodigal
	os.chdir(path_prodigal) 
	cmd="./prodigal.linux -i "+genome_complet+" -d "+cds_predits
	os.system(cmd)

	#2. Prediction 
	if scaffolds: 
		os.system(prediction_scaffold)
	if not scaffolds: 
		os.system(prediction) 

	#3. Blast 
	os.system(blast) 
			
		
	 
