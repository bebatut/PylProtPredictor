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
	        wx.Frame.__init__(self,parent,id,title)
	        self.parent = parent
		self.initialize() 
	
	#Design + actions#
	def initialize(self):
	        sizer = wx.GridBagSizer(5,15)
		
		espace_vertical_20=wx.Panel(self,-1,size=(0,20))
		espace_vertical_5=wx.Panel(self,-1,size=(0,5))

		#FONTS#
		titre=wx.Font(13,wx.NORMAL,wx.NORMAL,wx.BOLD)		
		sous_titre=wx.Font(10,wx.NORMAL,wx.NORMAL,wx.BOLD)
		font_valider=wx.Font(12,wx.NORMAL,wx.NORMAL,wx.BOLD)

		#DESIGN#
		sizer.Add(espace_vertical_20,(0,0),(1,5),wx.EXPAND)

		label_titre_blast=wx.StaticText(self,-1,label="BLAST")
		label_titre_blast.SetFont(titre)
		sizer.Add(label_titre_blast,(1,1),(1,5),wx.EXPAND)	

		line=wx.StaticLine(self,-1)
		sizer.Add(line,(2,0),(1,5),wx.EXPAND)

		sizer.Add(espace_vertical_5,(3,0),(1,5),wx.EXPAND)

		label_titre_exe=wx.StaticText(self,-1,label="Execution") 
		label_titre_exe.SetFont(sous_titre)
		sizer.Add(label_titre_exe,(4,1),(1,4),wx.EXPAND)	

		sizer.Add(espace_vertical_5,(5,0),(1,5),wx.EXPAND)

		label_prot_plus=wx.StaticText(self,-1,label="Potential PYL proteins (strand+) (fasta)" )
		sizer.Add(label_prot_plus,(6,1),(1,1),wx.EXPAND)
		filepicker_prot_plus=wx.FilePickerCtrl(self,-1) 
		sizer.Add(filepicker_prot_plus,(7,1),(1,1),wx.EXPAND)	

		label_prot_moins=wx.StaticText(self,-1,label="Potential PYL proteins (strand-) (fasta)")
		sizer.Add(label_prot_moins,(6,3),(1,1),wx.EXPAND) 
		filepicker_prot_moins=wx.FilePickerCtrl(self,-1)
		sizer.Add(filepicker_prot_moins,(7,3),(1,1),wx.EXPAND)

		sizer.Add(espace_vertical_5,(8,0),(1,5),wx.EXPAND)

		label_db=wx.StaticText(self,-1,label="Database") 
		sizer.Add(label_db,(9,1),(1,1),wx.EXPAND) 
		filepicker_db=wx.FilePickerCtrl(self,-1) 
		sizer.Add(filepicker_db,(10,1),(1,1),wx.EXPAND)
		
		label_evalue=wx.StaticText(self,-1,label="E-value max (if none, it will be 10)") 
		sizer.Add(label_evalue,(9,3),(1,1),wx.EXPAND) 
		text_evalue=wx.TextCtrl(self,-1)
		sizer.Add(text_evalue,(10,3),(1,1),wx.EXPAND)

		sizer.Add(espace_vertical_5,(11,0),(1,5),wx.EXPAND)

		checkbox_format=wx.CheckBox(self,-1,label="Database already formatted")
		sizer.Add(checkbox_format,(12,1),(1,1),wx.EXPAND)

		sizer.Add(espace_vertical_20,(13,0),(1,5),wx.EXPAND)

		titre_tri=wx.StaticText(self,-1,label="Out results") 
		titre_tri.SetFont(sous_titre) 
		sizer.Add(titre_tri,(14,1),(1,4),wx.EXPAND)
		
		sizer.Add(espace_vertical_5,(15,0),(1,5),wx.EXPAND)

		checkbox_align=wx.CheckBox(self,-1,label="Create simple file with alignment titles\nand only sequences for first alignment")
		sizer.Add(checkbox_align,(16,1),(1,1),wx.EXPAND)

		sizer.Add(espace_vertical_5,(17,0),(1,5),wx.EXPAND) 

		label_result=wx.StaticText(self,-1,label="Blast results directory") 
		sizer.Add(label_result,(18,1),(1,1),wx.EXPAND) 
		dirpicker_result=wx.DirPickerCtrl(self,-1)
		sizer.Add(dirpicker_result,(19,1),(1,1),wx.EXPAND)

		sizer.Add(espace_vertical_20,(20,0),(1,5),wx.EXPAND)

		button_ok=wx.Button(self,-1,label="OK") 
		button_ok.SetBackgroundColour("lightgrey") 
		button_ok.SetFont(font_valider) 
		sizer.Add(button_ok,(21,1),(1,3),wx.EXPAND)

		sizer.Add(espace_vertical_20,(22,0),(1,5),wx.EXPAND)

		self.Bind(wx.EVT_BUTTON,lambda event:self.verify_file(event,filepicker_prot_plus.GetPath(),filepicker_prot_moins.GetPath(),filepicker_db.GetPath(),dirpicker_result.GetPath(),checkbox_align.GetValue(),text_evalue.GetValue(),checkbox_format.GetValue()),button_ok) #Action : clic sur bouton OK = lancement launch_script, passe tous les arguments 

		self.SetSizerAndFit(sizer)
		self.Show(True)


	def verify_file(self,event,path_prot_plus,path_prot_moins,path_db,path_result,align,evalue,format):
		os.system("clear") 		
		prot_plus=True
		prot_moins=True
		db=True
		result=True
		bool_evalue=True
		
		if path_prot_plus=="": #pas de fichier 
			print "Prot plus : select a file" 
			prot_plus=False
		if path_prot_moins=="": 
			print "Prot moins : select a file" 
			prot_moins=False
		if path_db=="": 
			print "Database : select a file"
			db=False
		if path_result=="":
			print "Results directory : select a directory" 
			result=False
		if evalue=="": #si pas de valeur ds evalue 
			evalue="10" #valeur par defaut 
		expr="^\d+|^\d+e\d" #expression reguliere : x ou xex (où x=1 ou plusieurs chiffres) 
		if not (re.search(expr,evalue)): #si evalue ne correspond pas a l'expression reguliere 
			print "Evalue : need a number"  
			evalue_bool=False
		if prot_plus and prot_moins and db and result and evalue: #si tous les parametres OK 
			self.launch_script(event,path_prot_plus,path_prot_moins,path_db,path_result,align,evalue,format)
	       
	def launch_script(self,event,path_prot_plus,path_prot_moins,path_db,path_result,align,evalue,format):
		if not format: #formatage bdd 
			cmd="formatdb -i "+path_db
			print "format db..."
			os.system(cmd)
		home_directory=os.environ['HOME'] 
		script=home_directory+"/pyl_protein_prediction/Scripts/blast.py"  
		cmd="python "+script+" "+path_prot_plus+" "+path_prot_moins+" "+path_db+" "+path_result+" "+str(evalue)+" "+str(align)
		os.system(cmd)
		

#Main creation fenetre
if __name__ == "__main__":
    app = wx.App()  
    frame = simpleapp_wx(None,-1,'Pipeline:Blast') 
    app.MainLoop()
