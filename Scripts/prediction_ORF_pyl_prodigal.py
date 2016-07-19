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
				out_position_complement.write(position+"\n")  #ecrit la position dans le bon fichier 
			elif (strand=='1'): #sinon si brin plus 
				liste_seq_normal.append(record) #ajout cds a la bonne liste  
				out_position_normal.write(position+"\n") #ecrit la position dans le bon fichier 
	SeqIO.write(liste_seq_complement,sortie_cds_complement,"fasta") #ecrit la liste des CDS brin - dans un fichier fasta 
	SeqIO.write(liste_seq_normal,sortie_cds_normal,"fasta") #ecrit la liste des CDS brin + dans un fichier fasta  	

#Calcul taille genome 
def calculer_taille_genome(genome): 
	record=SeqIO.read(genome,"fasta") 
	return len(record)  

			
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

#Traduction (avc pyrrolysine) de tous les fichier fasta d'un dossier 
def traduction_tout_dossier(dossier,dossier_sortie_prot,dossier_sortie_prot_pyl): 
	for file in os.listdir(dossier):
		path_file=dossier+file
		if(os.path.isfile(path_file)):
			nom=file.replace(".fasta","") 
			path_sortie=dossier_sortie_prot+nom+"_protein.fasta"
			traduction_proteine(path_file,path_sortie) 
			trad_pyrrolysine(path_sortie,dossier_sortie_prot_pyl+nom+"_proteinPYL.fasta") 

#Recupere position allant du debut du debut du CDS jusqu'au debut du CDS suivant (brin +)
def position_region_a_conserver_normal(position_cds,fichier_sortie,taille_genome): 
	in_cds=open(position_cds,"r") #ouverture fichier position CDS 
	list_start=[] #initialisation liste début des CDS 
	list_end=[] #initialisation liste fin des CDS 
	out=open(fichier_sortie,"w") #ouverture en écriture fichier de sortie
	for position in in_cds: #parcours des positions du fichier des positions CDS 
		split_position=position.split("..") #split par rapport à .. pour récupérer start et end des CDS 
		start=split_position[0] 
		end=split_position[1] 
		list_start.append(start)  #ajout à la liste des start 
		list_end.append(end) #ajout à la liste des end 
	i=0 
	while i<len(list_start): #parcours des listes en même temps 
		start_region=list_start[i] #debut de la region = debut du CDS 
		if not i==len(list_start)-1: #si i n'est pas le dernier élément 
			if int(list_end[i])>int(list_start[i+1]): #si fin du CDS est plus grand que le debut du CDS suivant = CDS chevauchant 
				 end_region=int(list_start[i+2])-1 #fin de la region = debut du CDS +2 après (on saute le CDS chevauchant qui est directement en aval) -1 pour ne pas prendre le premier nt du codon start 
			else:
				end_region=int(list_start[i+1])-1 #si pas CDS chevauchant, fin de la region = début du CDS suivant (-1 pour ne pas prendre le premier nt du codon start)   
			position=start_region+".."+str(end_region) #ecriture propre position	
		elif i==len(list_start)-1: #si on est au dernier element du tableau
			end_region=taille_genome #fin de la region = fin du génome 
			position=start_region+".."+str(end_region)
		out.write(position+"\n") #ecriture des positions dans le fichier de sortie	 
		i+=1

#Recupere position allant de la fin du CDS a la fin du CDS suivant (brin -) 
def position_region_a_conserver_complement(position_cds,fichier_sortie):
	in_cds=open(position_cds,"r") 
	list_start=[]
	list_end=[]
	out=open(fichier_sortie,"w") 
	for position in in_cds:
		split_position=position.split("..") 
		start=split_position[0] 
		end=split_position[1] 
		list_start.append(start) 
		list_end.append(end)
	i=0 
	while i<len(list_start):
		end_region=list_end[i] #fin de la region = fin du CDS 
		if i!=0: #si i n'est pas le premier élément 
			if int(list_start[i])<int(list_end[i-1]): # si debut du CDS inferieur à la fin du CDS d'avant = CDS chevauchants  
				 start_region=int(list_end[i-2])+1 #debut de la region = fin de la region du CDS -2 en amont (on saute le CDS chevauchant) 
			else: #pas de CDS chevauchant 
				start_region=int(list_end[i-1])+1 #debut de la region = fin du CDS d'avant (+1 pour ne pas prendre le dernier nt du codon stop) 
			position=str(start_region)+".."+end_region	
		elif i==0: #si premiere iteration 
			start_region=1 #debut de la region = debut de la sequence genomique 
			position=str(start_region)+".."+end_region 		
		out.write(position)	 
		i+=1

#Recupere les regions a conserver à partir du génome complet
def recup_region_a_conserver(genome,fichier_position,fichier_sortie):
	record=SeqIO.read(genome,"fasta") #lecture sequence genomique  
	input_intergenique=open(fichier_position,"r") #lecture fichier positions 
	liste_region_a_conserver=[] #initialisation liste des regions 
	count=0 
	for position in input_intergenique: 
		count+=1 #incrementation a chaque nouvelle position 
		split_position=position.split("..") #split des positions pour separer debut/fin 
		start=split_position[0]
		end=split_position[1]
		end=end.replace("\n","") 
		region=record[int(start)-1:int(end)] #recuperation de la sequence sur le genome 
		region.description="CDS+intergenique "+str(count) #ajout a la description de "CDS+intergenique numero position"  
		liste_region_a_conserver.append(region) #ajout de l'objet SeqRecord a la liste   
	SeqIO.write(liste_region_a_conserver,fichier_sortie,"fasta") #ecriture du fasta de sortie 


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
	fasta1=False
	fasta2=False
	input_file=open(file,"r") 
	lire_file=input_file.read() 
	if (re.search(r"^>",lire_file)): 
		fasta1=True
	for record in SeqIO.parse(file,"fasta"):
		fasta2=True	
	if fasta1 and fasta2: 
		fasta=True 	  
	return fasta

#Verifie si le fichier contient plusieurs scaffolds 
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
CDS_TAG_plus=dir_result+"/CDS_TAG_plus.fasta" #Contient tous les CDS se terminant par TAG du brin+  
CDS_TAG_moins=dir_result+"/CDS_TAG_moins.fasta" #Contient tous les CDS se terminant par TAG du brin-
position_CDS_TAG_plus=dir_result+"/position_CDS_TAG_plus.txt" #Contient la position sur le genome des CDS TAG du brin+
position_CDS_TAG_moins=dir_result+"/position_CDS_TAG_moins.txt" #Contient la position sur le genome des CDS TAG du brin- 
position_cds_intergen_plus=dir_result+"/position_keep_regions_normal.txt" #Contient les positions des regions a conserver (CDS+region intergenique en aval) (brin+) 
position_cds_intergen_moins=dir_result+"/position_keep_region_complement.txt" #Contient les positions des regions a conserver (CDS+region intergenique en amont) (brin-) 
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

#Verification format fichier 

print "\nPREDICTION POTENTIAL PYL PROTEINS" 
print "_________________________\n"

print "Check files..." 
fasta_genome=verif_fasta(genome_complet) 
fasta_CDS=verif_fasta(CDS_predits) 
scaffold=verif_scaffold(genome_complet) 

if not fasta_genome: 
	print "Complete genome : need fasta format" 
	
if not fasta_CDS: 
	print "Predicted CDS : need fasta format"

if scaffold: 
	print "Complete genome : several scaffolds, select \"scaffolds\""
 	 
if fasta_genome and fasta_CDS and not scaffold:
	#2. Recuperer les cds se terminant par tag dans deux fichiers de sortie (un pour brin+, un pour brin-). Pareil pour les position
	if not os.path.exists(CDS_TAG_plus):
		print "Create fasta file with CDS TAG and position file of CDS TAG..."
		recuperer_cds_position_tag_tri(CDS_predits,CDS_TAG_plus,CDS_TAG_moins,position_CDS_TAG_plus,position_CDS_TAG_moins)

	#3. Recuperer les positions des regions a conserver (CDS+region intergenique en amont/aval selon le brin), pour ensuite les extraire du genome complet 
	if not os.path.exists(position_cds_intergen_plus):  
		print "Position regions to keep (+)..." 
		taille_genome=calculer_taille_genome(genome_complet)
		position_region_a_conserver_normal(position_CDS_TAG_plus,position_cds_intergen_plus,taille_genome)
	if not os.path.exists(position_cds_intergen_moins): 
		print "Position region to keep (-)..."
		position_region_a_conserver_complement(position_CDS_TAG_moins,position_cds_intergen_moins) 

	#4. Recuperer les sequences à conserver a partir du genome complet. Reverse-complementation des regions du brin -  
	if not os.path.exists(region_cds_intergen_plus): 
		print"Regions to keep (+)..."
		recup_region_a_conserver(genome_complet,position_cds_intergen_plus,region_cds_intergen_plus)
	if not os.path.exists(region_cds_intergen_moins):
		print"Regions to keep (-)..."  		
		recup_region_a_conserver(genome_complet,position_cds_intergen_moins,region_cds_intergen_moins)
		reverse_complementation_fichier_fasta(region_cds_intergen_moins,region_cds_intergen_moins_RC) 

	#5. Liste des genes possibles pour chaque CDS (avec ou sans pyl) 
	if not os.path.exists(genes_possibles_plus): 	
		print "Potential genes (+)..." 
		determiner_genes_possibles(region_cds_intergen_plus,genes_possibles_plus)
	if not os.path.exists(genes_possibles_moins): 
		print "Potential genes (-)..." 
		determiner_genes_possibles(region_cds_intergen_moins_RC,genes_possibles_moins)

	#6. Traduction des genes (avec et sans pyl) 
	if not os.path.exists(protein_no_pyl_plus):
		print "Traduction (+)..." 
		traduction_proteine(genes_possibles_plus,protein_no_pyl_plus)
		trad_pyrrolysine(protein_no_pyl_plus,protein_pyl_plus)  
	if not os.path.exists(protein_no_pyl_moins): 
		print "Traduction (-)..."
		traduction_proteine(genes_possibles_moins,protein_no_pyl_moins) 
		trad_pyrrolysine(protein_no_pyl_moins,protein_pyl_moins)

	print "Prediction done!"	
