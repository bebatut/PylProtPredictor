#!/usr/bin/python2.7
# -*- coding: utf-8 -*- 

from Bio import SeqIO 
from Bio import SeqRecord 
from Bio.Seq import Seq 
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_protein
import re 
import os 
import sys 

##########
#METHODES#
##########

#Dans le fichier .fasta contenant tous les CDS predits, recupere uniquement les CDS se terminant par TAG et les trie en fonction du brin sur lequel ils se trouvent 
def recuperer_cds_position_tag_tri(fichier_cds,sortie_cds_normal,sortie_cds_complement,sortie_position_normal,sortie_position_complement):
	out_position_normal=open(sortie_position_normal,"w") 
	out_position_complement=open(sortie_position_complement,"w") 	
	liste_seq_complement=[] #initialise liste cds du brin -
	liste_seq_normal=[] #initialise liste cds du brin + 
	for record in SeqIO.parse(fichier_cds,"fasta"): #parcours du fichier de tous les CDS 
		if (re.search(r"TAG$",str(record.seq))): #si se termine par TAG 
			split_description=record.description.split("#") #split de la description pour recuperer les positions 
			start=split_description[1] #recuperation position depart 
			start=start.replace(" ","") #enleve les espaces 
			end=split_description[2] #recuperation position fin 
			end=end.replace(" ","") #enleve les espace 
			strand=split_description[3] #recupere le brin (1 ou -1) 
			strand=strand.replace(" ","") #enleve espaces 
			position=start+".."+end #position entiere forme : start..end 
			if (strand=='-1'): #si brin moins 
				liste_seq_complement.append(record) # ajoute le cds a la bonne liste 
				out_position_complement.write(record.id+" "+position+"\n")  #ecrit la position dans le bon fichier 
			elif (strand=='1'): #sinon si brin plus 
				liste_seq_normal.append(record) #ajout cds a la bonne liste  
				out_position_normal.write(record.id+" "+position+"\n") #ecrit la position dans le bon fichier 
	SeqIO.write(liste_seq_complement,sortie_cds_complement,"fasta") #ecrit la liste des CDS brin - dans un fichier fasta 
	SeqIO.write(liste_seq_normal,sortie_cds_normal,"fasta") #ecrit la liste des CDS brin + dans un fichier fasta  	

#Calcul taille scaffolds et les ecrire dans un fichier de sortie
def calculer_taille_scaffolds(genome,out_file): 
	out_file=open(out_file,"w") 
	for record in SeqIO.parse(genome,"fasta"): 
		out_file.write(record.id+" "+str(len(record))+"\n")  
	  

			
#Creation d'une liste des codons d'une sequence (separe la seq par 3) 
def creer_liste_codons(sequence):
	list_codon=[]
	while (len(sequence)!=0): #tant que la sequence contient quelque chose 
		codon=sequence[:3] #recupere les 3 premiers nt 
		sequence=sequence[3:] #elimine les 3 premiers nt 
		if (len(codon)==3): #ajoute a la liste uniquement si le codon contient 3 nt 
			codon_str=str(codon.seq)
			list_codon.append(codon_str)
	return list_codon

#Recherche du premier codon stop a partir d'une liste de codon 
def recherche_codon_stop(list_codon):
	continuer='vrai' 
	i=0
	while(continuer=='vrai'): 
		while i<len(list_codon):
			if list_codon[i]=='TAG': 
				continuer='faux' 
				return i,"tag" 
			elif list_codon[i]=='TAA' or list_codon[i]=='TGA': 
				continuer='faux' 
				return i,"taa,tga" 
			i+=1 	
		return -1,"no_stop"

#Traduction en prot d'un fichier fasta nucleique 
def traduction_proteine(fichier_fasta,fichier_sortie):
	liste_prot=[] 
	for record in SeqIO.parse(fichier_fasta,"fasta"): 
		seq_prot=record.seq.translate() 
		record_prot=SeqRecord(seq_prot,id=record.id,description=record.description) 
		liste_prot.append(record_prot) 
	SeqIO.write(liste_prot,fichier_sortie,"fasta")

#Traduction en prot en prenant en compte que les TAG pas a la fin codent pour pyrrolysine 
def trad_pyrrolysine(fasta_prot,sortie): 
	liste_prot=[]	
	for record in SeqIO.parse(fasta_prot,"fasta"): 
		prot=''
		prot=str(record.seq).replace("*","O") 
		prot=prot[:-1]
		prot+="*" 
		prot_seq=Seq(prot,generic_protein) 
		record_prot=SeqRecord(prot_seq,id=record.id,description=record.description) 
		liste_prot.append(record_prot) 
	SeqIO.write(liste_prot,sortie,"fasta") 

#Recupère position allant du début du CDS au début du CDS suivant 
#Pour dernier CDS d'un scaffold, fin = taille du scaffold
def position_region_a_conserver_normal(position_cds,fichier_sortie,taille_scaffolds): 
	in_cds=open(position_cds,"r") #ouverture fichier position CDS 
	dict_start_cds={} #initialisation tableau association scaffold ->debut cds  
	dict_end_cds={} #initialisation tableau association scaffold->fin cds 
	dict_taille={} #initialisation tableau association scaffold->taille scaffold
	dict_number_scaffold={}
	taille_scaffolds=open(taille_scaffolds,"r") #ouverture fichier taille scaffold 
	out=open(fichier_sortie,"w") #ouverture en écriture fichier de sortie

	for line in taille_scaffolds: 
		split_line=line.split(" ") 
		taille=split_line[1]
		taille=taille.replace("\n","") 
		scaffold=split_line[0]
		dict_taille[scaffold]=taille 
		dict_start_cds[scaffold]=[]
		dict_end_cds[scaffold]=[]
		dict_number_scaffold[scaffold]=[]

	for line in in_cds: 
		split_line=line.split(" ") 
		position=split_line[1]
		split_position=position.split("..") 
		start=split_position[0]
		end=split_position[1]
		position=position.replace("\n","") 
		scaffold=split_line[0]
		split_scaffold=scaffold.split("_") 
		number_scaffold=split_scaffold[2]
		scaffold=split_scaffold[0]+"_"+split_scaffold[1]
		dict_number_scaffold[scaffold].append(number_scaffold) 
		dict_start_cds[scaffold].append(start)
		dict_end_cds[scaffold].append(end)

	for scaffold, start in dict_start_cds.items(): 
		if not len(dict_start_cds[scaffold])==0: 
			for i,start in enumerate(dict_start_cds[scaffold]): #parcours de toutes les valeurs du dico pr chaque clé
				region_start=dict_start_cds[scaffold][i] 
				if not i == len(dict_start_cds[scaffold])-1: #si on est pas à la derniere iteration
					if int(dict_end_cds[scaffold][i])>int(dict_start_cds[scaffold][i+1]): #si fin du CDS > debut du CDS d'apres = CDS chevauchant
						if i==len(dict_start_cds[scaffold])-2: #si il y a seulement un CDS en amont
							region_end=dict_taille[scaffold] #fin de la region=fin scafold
						else: #s'il y a plus d'1 cds en amont, on peut donc aller chercher le CDS à la position +2
							region_end=int(dict_start_cds[scaffold][i+2])-1 #fin de la region=debut CDS +2 en aval 
					else: #si pas chevauchant 
						region_end=int(dict_start_cds[scaffold][i+1])-1 #fin de la region=debut CDS en aval 
				else: #si on est à la derniere itération 
					region_end=dict_taille[scaffold] #fin=fin scaffold(sa taille)
				region_end=str(region_end)
	 			position=region_start+".."+region_end
				scaffold_complet=scaffold+"_"+dict_number_scaffold[scaffold][i]
				out.write(scaffold_complet+" "+position+"\n")


#Recupere position allant de la fin du CDS a la fin du CDS suivant (brin -)
#Pour premier CDS d'un scaffold, debut=1 
def position_region_a_conserver_complement(position_cds,fichier_sortie,taille_scaffolds):
	in_cds=open(position_cds,"r") 
	dict_end_cds={} #dictionnaire scaffold (scaffold_X) -> fin_cds
	dict_start_cds={} #dictionnaire scaffold -> debut_cds 
	dict_number_scaffold={} #dictionnaire scaffold -> numeros des differents CDS (x dans scaffold_X_x) 
	out=open(fichier_sortie,"w") 
	taille_scaffolds=open(taille_scaffolds,"r") 
		
	for line in taille_scaffolds:  
		split_line=line.split(" ") 
		scaffold=split_line[0] 
		dict_end_cds[scaffold]=[] 
		dict_start_cds[scaffold]=[]
		dict_number_scaffold[scaffold]=[]		

	for line in in_cds: #parcours des positions du fichier des positions CDS 
		split_line=line.split(" ") #split par rapport à .. pour récupérer start et end des CDS 
		scaffold=split_line[0]
		split_scaffold=scaffold.split("_") 
		number_scaffold=split_scaffold[2]
		scaffold=split_scaffold[0]+"_"+split_scaffold[1]
		position=split_line[1]
		position=position.split("..")
		start=position[0] 
		end=position[1] 
		end=end.replace("\n","") 
		dict_end_cds[scaffold].append(end)
		dict_start_cds[scaffold].append(start)
		dict_number_scaffold[scaffold].append(number_scaffold)

	for scaffold,values in dict_end_cds.items(): 
		if not len(dict_end_cds[scaffold])==0: #si le scaffold a des valeurs 
			for i,values in enumerate(dict_end_cds[scaffold]): 
				region_end=dict_end_cds[scaffold][i] #fin de la region=fin du CDS
				if not i==0: #si on est pas a la 1e iteration 
					if int(dict_start_cds[scaffold][i])<int(dict_end_cds[scaffold][i-1]):#si debut du CDS < fin du CDS d'avant = chevauchant 
						if i==1: #si il y a seulement un CDS en amont
							region_start=1 #debut_region=1
						else: 
							region_start=int(dict_end_cds[scaffold][i-2])+1 #debut region=fin CDS -2 en amont 
					else: 
						region_start=int(dict_end_cds[scaffold][i-1])+1 #debut region=fin CDS en amont 
					
				else: #si on est a la 1ere iteration
					region_start=1 #debut region=1
				region_start=str(region_start)
				position=region_start+".."+region_end
				scaffold_complet=scaffold+"_"+dict_number_scaffold[scaffold][i]
				out.write(scaffold_complet+" "+position+"\n")
				
	

#Recupere les regions a conserver à partir du génome complet
def recup_region_a_conserver(genome,fichier_position,fichier_sortie):
	liste_region_a_conserver=[] 
	position_a_conserver=open(fichier_position,"r")
	for line in position_a_conserver: 
		split_line=line.split(" ") 
		scaffold_entier=split_line[0] #scaffold_X_x
		scaffold=scaffold_entier.split("_") 
		scaffold=scaffold[0]+"_"+scaffold[1] #scaffold_X
		position=split_line[1]
		split_position=position.split("..") 
		start=split_position[0]
		end=split_position[1]
		end=end.replace("\n","")
		pattern=scaffold+"$" #pour expression reguliere, recherche id se terminant par scaffold_X (pour ne pas avoir par ex scaffold_Xxxx comme resultat)
		for record in SeqIO.parse(genome,"fasta"):
			if (re.search(pattern,record.id)):
				region=record[int(start)-1:int(end)] #recuperation region 
				region.id=scaffold_entier
				region.description=""
				liste_region_a_conserver.append(region) 
		SeqIO.write(liste_region_a_conserver,fichier_sortie,"fasta") 
	

#Reverse complemente un fichier fasta entier 
def reverse_complementation_fichier_fasta(fichier_fasta,fichier_sortie): 
	liste_record_rc=[]
	for record in SeqIO.parse(fichier_fasta,"fasta"): 
		record_rc=record.reverse_complement()
		record_rc.id=record.id 
		record_rc.description=record.description
		liste_record_rc.append(record_rc)  
	SeqIO.write(liste_record_rc,fichier_sortie,"fasta") 

#Determine tous les genes possibles (codon tag comme codon stop ou codant, si 2e "stop" = tag, "sauvegarde" du gene et prolongation jusqu'au prochain codon stop
def determiner_genes_possibles(region_intergen_cds,fichier_sortie):
	liste_gene=[]
	i=0
	for record in SeqIO.parse(region_intergen_cds,"fasta"): #parcours fichier region conservee
		i+=1 
		count_tag=0 #nombre tag 			
		liste_codon=creer_liste_codons(record) #separe la seq en codon 
		position_codon_stop,codon=recherche_codon_stop(liste_codon) #recherche position codon stop et nature
		if codon=="no_stop": #si on ne trouve pas de codon stop 
			continue #prochaine iteration de for 
		position_dernier_nt_stop=(position_codon_stop+1)*3 #donne position du dernier nt prcq on a juste la position du codon 
		gene=record[:position_dernier_nt_stop] #recuperation de la seq du debut au dernier nt du codon stop 
		gene.description=record.description.replace("+intergenique","")+" NO PYL" 
		liste_gene.append(gene) #remplissage tableau genes 
		tag='true' 
		while(tag=='true'): #commence avec le 1er codon qui est forcement un tag  
			count_tag+=1 #incrementation compteur tag 
			liste_codon=liste_codon[position_codon_stop+1:] #coupage du tableau des codons (seq) du premier codon stop a la fin 
			position_codon_stop,codon=recherche_codon_stop(liste_codon) #recherche du prochain codon dans le nouveau tableau
			if codon=="no_stop": #si on ne trouve pas de codon stop 
				break #sortie de la boucle while
			position_dernier_nt_stop_2=(position_codon_stop+1)*3 
			position_dernier_nt_stop+=position_dernier_nt_stop_2 #position dernier nt = position dernier nt dans le tableau + position dernier nt d'avant 
			gene=record[:position_dernier_nt_stop] #decoupage du gene du debut au dernier nt 
			gene.description=record.description.replace("+intergenique","")+" PYL "+str(count_tag)
			liste_gene.append(gene) #ajout gene  
			if codon=="taa,tga": #si codon est taa ou tga 
				tag='false' 
	SeqIO.write(liste_gene,fichier_sortie,"fasta") #ecriture de tous les genes dans fichier fasta   

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



			 

	
######
#MAIN#
######


#GESTION DES FICHIERS#
genome_complet=sys.argv[1]
CDS_predits=sys.argv[2]
dir_result=sys.argv[3]
CDS_TAG_plus=dir_result+"/CDS_TAG_normal.fasta" #Contient tous les CDS se terminant par TAG du brin+  
CDS_TAG_moins=dir_result+"/CDS_TAG_complement.fasta" #Contient tous les CDS se terminant par TAG du brin-
position_CDS_TAG_plus=dir_result+"/position_CDS_TAG_normal.txt" #Contient la position sur le genome des CDS TAG du brin+
position_CDS_TAG_moins=dir_result+"/position_CDS_TAG_complement.txt" #Contient la position sur le genome des CDS TAG du brin- 
fichier_longueur_scaffolds=dir_result+"/length_scaffolds.txt" 
position_cds_intergen_plus=dir_result+"/position_keep_regions_normal.txt" #Contient les positions des regions a conserver (CDS+region intergenique en aval) (brin+) 
position_cds_intergen_moins=dir_result+"/position_keep_regions_complement.txt" #Contient les positions des regions a conserver (CDS+region intergenique en amont) (brin-) 
region_cds_intergen_plus=dir_result+"/keep_regions_normal.fasta" #Contient les regions a conserver (brin+)  
region_cds_intergen_moins=dir_result+"/keep_regions_complement.fasta" #Contient les regions a conserver (brin-) 
region_cds_intergen_moins_RC=dir_result+"/keep_regions_complement_RC.fasta" #Regions reverse complementees pour le brin -
genes_possibles_plus=dir_result+"/potential_genes_normal.fasta" #Liste des genes possibles (brin+) 
genes_possibles_moins=dir_result+"/potential_genes_complement.fasta" #Liste des genes possibles (brin-) 
protein_no_pyl_plus=dir_result+"/potential_proteins_normal.fasta" #Genes possibles traduits (brin+) 
protein_no_pyl_moins=dir_result+"/potential_proteins_complement.fasta"#Genes possibles traduits (brin-) 
protein_pyl_plus=dir_result+"/potential_PYL_proteins_normal.fasta" #Genes possibles traduits avec les tag codants = O (brin+) 
protein_pyl_moins=dir_result+"/potential_PYL_proteins_complement.fasta" #Genes possibles traduits avec les tag codants = O (brin-)


#EXECUTION#


print "\nPREDICTION POTENTIAL PYL PROTEINS" 
print "_________________________\n"

#Verification fichier entree 
print "Check files..." 
fasta_genome=verif_fasta(genome_complet) 
fasta_CDS=verif_fasta(CDS_predits) 
scaffolds=verif_scaffold(genome_complet)
if not fasta_genome: 
	print "Complete genome : need fasta format" 
if not fasta_CDS: 
	print "Predicted CDS : need fasta format" 
if not scaffolds : 
	print "Complete genome : genome is full assembled, select \"fully assembled genome\""
 
if fasta_genome and fasta_CDS and scaffolds: 
	#1. Recuperer les cds se terminant par tag dans deux fichiers de sortie (un pour brin+, un pour brin-). Pareil pour les position 
	if not os.path.exists(CDS_TAG_plus):
		print "Create fasta file with CDS TAG and position file of CDS TAG..."
		recuperer_cds_position_tag_tri(CDS_predits,CDS_TAG_plus,CDS_TAG_moins,position_CDS_TAG_plus,position_CDS_TAG_moins)

	#2. Recuperer les positions des regions a conserver (CDS+region intergenique en amont/aval selon le brin), pour ensuite les extraire du genome complet 
	calculer_taille_scaffolds(genome_complet,fichier_longueur_scaffolds) 
	if not os.path.exists(position_cds_intergen_plus):  
		print "Position regions to keep (+)..." 
		position_region_a_conserver_normal(position_CDS_TAG_plus,position_cds_intergen_plus,fichier_longueur_scaffolds)
	if not os.path.exists(position_cds_intergen_moins): 
		print "Position regions to keep (-)..."
		position_region_a_conserver_complement(position_CDS_TAG_moins,position_cds_intergen_moins,fichier_longueur_scaffolds) 

	#3. Recuperer les sequences à conserver a partir du genome complet. Reverse-complementation des regions du brin -  
	if not os.path.exists(region_cds_intergen_plus): 
		print"Regions to keep (+)..."
		recup_region_a_conserver(genome_complet,position_cds_intergen_plus,region_cds_intergen_plus)
	if not os.path.exists(region_cds_intergen_moins):
		print"Regions to keep (-)..."  		
		recup_region_a_conserver(genome_complet,position_cds_intergen_moins,region_cds_intergen_moins)
		reverse_complementation_fichier_fasta(region_cds_intergen_moins,region_cds_intergen_moins_RC) 

	#4. Liste des genes possibles pour chaque CDS (avec ou sans pyl) 
	if not os.path.exists(genes_possibles_plus): 	
		print "Potential genes (+)..." 
		determiner_genes_possibles(region_cds_intergen_plus,genes_possibles_plus)
	if not os.path.exists(genes_possibles_moins): 
		print "Potential genes  (-)..." 
		determiner_genes_possibles(region_cds_intergen_moins_RC,genes_possibles_moins)

	#5. Traduction des genes (avec et sans pyl) 
	if not os.path.exists(protein_no_pyl_plus):
		print "Traduction (+)..." 
		traduction_proteine(genes_possibles_plus,protein_no_pyl_plus)
		trad_pyrrolysine(protein_no_pyl_plus,protein_pyl_plus)  
	if not os.path.exists(protein_no_pyl_moins): 
		print "Traduction (-)..."
		traduction_proteine(genes_possibles_moins,protein_no_pyl_moins) 
		trad_pyrrolysine(protein_no_pyl_moins,protein_pyl_moins)

	print "Prediction done!"	
