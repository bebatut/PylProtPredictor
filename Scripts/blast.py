#!/usr/bin/python2.7 

import os 
import sys 
from Bio.Blast import NCBIXML
import re
from Bio import SeqIO 

##########
#METHODES#
##########

def lancement_blast(fichier_fasta,fichier_sortie,type_blast,database,evalue): 
	blast_cmd="blastall -p "+type_blast+" -d "+database+" -i "+fichier_fasta+" -o "+fichier_sortie+" -e "+evalue+" -F F"  
	os.system(blast_cmd)

def lancement_blast_xml(fichier_fasta,fichier_sortie,type_blast,database,evalue): 
	blast_cmd="blastall -p "+type_blast+" -d "+database+" -i "+fichier_fasta+" -o "+fichier_sortie+" -m 7 -e "+evalue+" -F F" 
	os.system(blast_cmd)

#Recupere la liste des alignements pour chaque query et les sequences du meilleur alignement 
def recupere_prot(resultat_blast,fichier_sortie,resume_sortie):
	result_handle=open(resultat_blast)
	out=open(fichier_sortie,"w") 
	out_resume=open(resume_sortie,"w")   
	for blast_records in NCBIXML.parse(result_handle): 
		out.write("query="+blast_records.query+"\n"+"("+str(blast_records.query_letters)+" letters)\n\n")
		out.write("match protein\t\t\t\t\t\t\t\t\t\t\tscore\tevalue\n")
		for i,alignment in enumerate(blast_records.alignments): 
			for hsp in alignment.hsps: 
				out.write(alignment.hit_def+"\t"+str(hsp.score)+"\t"+str(hsp.expect)+"\n") #liste des alignements avec def score evalue  
			if i==0: #seulement le premier alignement 
				for hsp in alignment.hsps:
					out_resume.write(blast_records.query+"\t"+str(hsp.align_length)+"\n")
					out.write("length="+str(alignment.length)+"\n") 
					alignement=alignment_format(hsp.query,hsp.match,hsp.sbjct) #retourne sequence alignement 
					out.write("\n"+alignement)
		out.write("_________________________________________________\n\n")


#Verifie si le fichier est sous format fasta 							
def verif_fasta(file): 
	fasta=False
	for record in SeqIO.parse(file,"fasta"):
		fasta=True		  
	return fasta

#Formate l'alignement pour affichage comme une sortie blast
def alignment_format(query,match,sbjct): 
	alignement=""
	while query !="" and match !="" and sbjct !="": 
		query_line=query[:60]
		match_line=match[:60]
		sbjct_line=sbjct[:60]
		alignement=alignement+query_line+"\n"+match_line+"\n"+sbjct_line+"\n\n"
		query=query[60:]
		match=match[60:]
		sbjct=sbjct[60:]
	return alignement

######					
#MAIN#
######

#GESTION DES FICHIERS#
dossier_result=sys.argv[4]
fichier_prot_normal=sys.argv[1]
fichier_prot_complement=sys.argv[2] 
evalue=sys.argv[5]
align=sys.argv[6]
resultat_blast_normal=dossier_result+"/blast_normal.out" #resultat "brut" blast brin +
resultat_blast_normal_xml=dossier_result+"/blast_normal.xml" #resultat "brut" xml (pr traitement BioPython) blast brin + 
resultat_blast_complement=dossier_result+"/blast_complement.out" #resultat brut blast brin - 
resultat_blast_complement_xml=dossier_result+"/blast_complemet.xml" #resultat brut xml (pr traitement BioPython) blast brin - 
traitement_blast_normal=dossier_result+"/simple_results_blast_normal.txt" #resultats simplifies brin + (liste alignements + sequence alignement seulement pour le 1er) 
traitement_blast_complement=dossier_result+"/simple_results_blast_complement.txt" #resultats simplifies brin -
resume_blast_normal=dossier_result+"/length_alignment_plus.txt" 
resume_blast_complement=dossier_result+"/length_alignment_moins.txt" 
db=sys.argv[3]
 

#EXECUTION#

print "\nBLASTP PYL PROTEINS"
print"__________________\n"

#Verification fichier entree 
print "Check files..." 
fasta_plus=verif_fasta(fichier_prot_normal)
fasta_moins=verif_fasta(fichier_prot_complement)
fasta_db=verif_fasta(db) 
if not fasta_plus: 
	print "Potentiel pyl proteins (strand+) : Need fasta format" 
if not fasta_moins: 
	print "Potentiel pyl proteins (strand-) : Need fasta format"
if not fasta_db: 
	print "Database : Need fasta format"  


if fasta_plus and fasta_moins:

	#1. Execution des blast 
	if not os.path.exists(resultat_blast_normal): 
		print "Classic blast (+)..." 
		lancement_blast(fichier_prot_normal,resultat_blast_normal,"blastp",db,evalue)
	if not os.path.exists(resultat_blast_complement): 
		print "Classic blast (-)..." 
		lancement_blast(fichier_prot_complement,resultat_blast_complement,"blastp",db,evalue)
	if not os.path.exists(resultat_blast_normal_xml):
		print "XML blast (+)..." 
		lancement_blast_xml(fichier_prot_normal,resultat_blast_normal_xml,"blastp",db,evalue) 	
	if not os.path.exists(resultat_blast_complement_xml):
		print "XML blast (-)..."
		lancement_blast_xml(fichier_prot_complement,resultat_blast_complement_xml,"blastp",db,evalue)

	#2.Traitement/resume des resultats
	if align=="True": 
		if not os.path.exists(resume_blast_normal):
			print "Results treatment (+)..." 
			recupere_prot(resultat_blast_normal_xml,traitement_blast_normal,resume_blast_normal)
		if not os.path.exists(resume_blast_complement): 
			print "Results treatment (-)..."  
			recupere_prot(resultat_blast_complement_xml,traitement_blast_complement,resume_blast_complement) 

	print "Blast done !" 
