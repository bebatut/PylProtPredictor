menu.py 
_____________________________________
Only this script will be "manually" launch 
Display a menu with choices :
- Prediction CDS prodigal : only CDS prediction with prodigal  
- Prediction pyl proteins : only pyl proteins prediction, need predicted CDS 
- Blast : only blast on predicted pyl proteins, need predicted proteins
- Complete pipeline : Prodigal+prediction pyl proteins+blast 


prodigal.py
_________________________________________
Only launch CDS prediction with prodigal 
Use the command : prodigal.linux with options -i (input file,complete genome) and -d (create a fasta file with CDS) 


prediction_pyl.py 
__________________________________________
Only launch pyl protein prediction 
Use the script : prediction_ORF_pyl_prodigal.py or prediction_ORF_pyl_prodigal_scaffold.py  


blast_pyl.py 
____________________________________________
Only launch blast and treatment on blast output 
Use the script : blast.py 


pipeline.py
____________________________________________
Launch all pipeline : prodigal+prediction pyl+blast 


