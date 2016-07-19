#!/usr/bin/python
# -*- coding: iso-8859-1 -*-
# -*- coding: utf-8 -*-
import sys 
import os 
import re 

try:
    import wx
except ImportError:
    raise ImportError,"The wxPython module is required to run this program"

class simpleapp_wx(wx.Frame):
	def __init__(self,parent,id,title):
	        wx.Frame.__init__(self,parent,id,title,size=(500,350))
	        self.parent = parent
		self.initialize()
	

	def initialize(self):
	        sizer = wx.GridBagSizer(5,15)

		espace_vertical_20=wx.Panel(self,-1,size=(0,20))
		espace_vertical_5=wx.Panel(self,-1,size=(0,5))
		
		#Fonts 
		titre=wx.Font(13,wx.NORMAL,wx.NORMAL,wx.BOLD)		
		sous_titre=wx.Font(10,wx.NORMAL,wx.NORMAL,wx.BOLD)
		font_valider=wx.Font(12,wx.NORMAL,wx.NORMAL,wx.BOLD)

		sizer.Add(espace_vertical_20,(0,0),(1,5),wx.EXPAND)

		label_titre_prediction=wx.StaticText(self,-1,label="PREDICTION")
		label_titre_prediction.SetFont(titre)
		sizer.Add(label_titre_prediction,(1,1),(1,4),wx.EXPAND)	

		label_ss_titre_prediction=wx.StaticText(self,-1,label="CDS/POTENTIAL PYL GENES/\nPOTENTIAL PYL PROTEINS") 
		label_ss_titre_prediction.SetFont(sous_titre)
		sizer.Add(label_ss_titre_prediction,(2,1),(1,4),wx.EXPAND)
		
		line=wx.StaticLine(self,-1)
		sizer.Add(line,(3,0),(1,5),wx.EXPAND)

		label_genome=wx.StaticText(self,-1,label="Complete genome (fasta)") 
		sizer.Add(label_genome,(5,1),(1,1),wx.EXPAND) 
		filepicker_genome=wx.FilePickerCtrl(self,-1) 
		sizer.Add(filepicker_genome,(6,1),(1,1),wx.EXPAND) 

		label_result_prediction=wx.StaticText(self,-1,label="Prediction results directory") 
		sizer.Add(label_result_prediction,(5,3),(1,1),wx.EXPAND)
		dirpicker_result_prediction=wx.DirPickerCtrl(self,-1) 
		sizer.Add(dirpicker_result_prediction,(6,3),(1,1),wx.EXPAND) 

		sizer.Add(espace_vertical_5,(7,0),(1,5),wx.EXPAND)
		line=wx.StaticLine(self,-1)
		sizer.Add(line,(8,0),(1,5),wx.EXPAND)
		line=wx.StaticLine(self,-1) 
		sizer.Add(line,(9,0),(1,5),wx.EXPAND)
		sizer.Add(espace_vertical_5,(10,0),(1,5),wx.EXPAND)

		label_titre_blast=wx.StaticText(self,-1,label="BLAST")
		label_titre_blast.SetFont(titre)
		sizer.Add(label_titre_blast,(11,1),(1,5),wx.EXPAND)	

		line=wx.StaticLine(self,-1)
		sizer.Add(line,(12,0),(1,5),wx.EXPAND)

		sizer.Add(espace_vertical_5,(13,0),(1,5),wx.EXPAND)

		label_db=wx.StaticText(self,-1,label="Database") 
		sizer.Add(label_db,(14,1),(1,1),wx.EXPAND) 
		filepicker_db=wx.FilePickerCtrl(self,-1) 
		sizer.Add(filepicker_db,(15,1),(1,1),wx.EXPAND)

		sizer.Add(espace_vertical_5,(16,0),(1,5),wx.EXPAND)

		checkbox_format=wx.CheckBox(self,-1,label="Database already formatted")
		sizer.Add(checkbox_format,(17,1),(1,1),wx.EXPAND)

		titre_tri=wx.StaticText(self,-1,label="Out results") 
		titre_tri.SetFont(sous_titre)		
		sizer.Add(titre_tri,(14,3),(1,1),wx.EXPAND)

		checkbox_align=wx.CheckBox(self,-1,label="Create simple file with alignment titles\nand only sequences for first alignment")
		sizer.Add(checkbox_align,(15,3),(1,1),wx.EXPAND)	
		
		sizer.Add(espace_vertical_5,(18,0),(1,1),wx.EXPAND)

		label_evalue=wx.StaticText(self,-1,label="E-value max (if none, it will be 10)") 
		sizer.Add(label_evalue,(19,1),(1,1),wx.EXPAND) 
		text_evalue=wx.TextCtrl(self,-1)
		sizer.Add(text_evalue,(20,1),(1,1),wx.EXPAND)

		sizer.Add(espace_vertical_5,(21,0),(1,5),wx.EXPAND) 

		label_result_blast=wx.StaticText(self,-1,label="Blast results directory") 
		sizer.Add(label_result_blast,(18,3),(1,1),wx.EXPAND) 
		dirpicker_result_blast=wx.DirPickerCtrl(self,-1)
		sizer.Add(dirpicker_result_blast,(19,3),(1,1),wx.EXPAND)

		sizer.Add(espace_vertical_5,(24,0),(1,5),wx.EXPAND)
		line=wx.StaticLine(self,-1)
		sizer.Add(line,(25,0),(1,5),wx.EXPAND)
		line=wx.StaticLine(self,-1) 
		sizer.Add(line,(26,0),(1,5),wx.EXPAND)
		sizer.Add(espace_vertical_5,(25,0),(1,5),wx.EXPAND)

		button_ok=wx.Button(self,-1,label="OK") 
		button_ok.SetBackgroundColour("lightgrey") 
		button_ok.SetFont(font_valider) 
		sizer.Add(button_ok,(27,1),(1,3),wx.EXPAND)

		sizer.Add(espace_vertical_20,(28,0),(1,5),wx.EXPAND)

		self.Bind(wx.EVT_BUTTON,lambda event: self.choice_script(event,filepicker_genome.GetPath(),dirpicker_result_prediction.GetPath(),filepicker_db.GetPath(),text_evalue.GetValue(),checkbox_align.GetValue(),dirpicker_result_blast.GetPath(),checkbox_format.GetValue()),button_ok) #lance choice_script qd clic sur bouton ok, passe tous les arguments 

		self.SetSizerAndFit(sizer)
		self.Show(True)

	def launch_pipeline(self,event,path_genome,path_result_prediction,path_db,evalue,align,path_result_blast,format):
		if not format: 
			cmd="formatdb -i "+path_db
			print "format db..."
			os.system(cmd)
		home_directory=os.environ['HOME']
		pipeline=home_directory+"/pyl_protein_prediction/Scripts/launch_pipeline.py "+path_genome+" "+path_result_prediction+" "+path_result_blast+" "+path_db+" "+evalue+" "+str(align)
		os.system(pipeline)
		self.Refresh(True)
		

	def choice_script(self,event,path_genome,path_result_prediction,path_db,evalue,align,path_result_blast,format):
		os.system("clear") 		
		genome=True
		results_prediction=True
		results_blast=True
		bool_evalue=True 
		
		if path_genome=="": 
			print "Complete genome : select a file" 
			genome=False
		if path_result_prediction=="": 
			print "Prediction results directory : select a directory"
			results_prediction=False
		if path_result_blast=="":
			print "Blast results directory : select a directory" 
			results_blast=False
		if evalue=="": 
			evalue="10" #valeur par defaut 
		expr="^\d+|^\d+e\d" #expression reguliere x ou xex (avec x=1 ou plusieurs chiffres) 
		if not (re.search(expr,evalue)):
			print "Evalue : need a number"  
			evalue_bool=False  
		if genome and results_prediction and results_blast and evalue: #si tous paramètres OK   
			self.launch_pipeline(event,path_genome,path_result_prediction,path_db,evalue,align,path_result_blast,format)
			self.Destroy() 
					
	
if __name__ == "__main__":
    app = wx.App()
    frame = simpleapp_wx(None,-1,'Pipeline:Complet')
    app.MainLoop()
