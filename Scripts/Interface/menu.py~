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
	        wx.Frame.__init__(self,parent,id,title,size=(400,170))
	        self.parent = parent
		self.initialize()
	
	#initialisation de la fenetre : design+actions sur les boutons 
	def initialize(self):
	        sizer = wx.GridBagSizer(20,10)

		panel=wx.Panel(self,-1)

		sizer.Add(panel,(1,0),(1,1),wx.EXPAND)

		button_prodigal=wx.Button(self,-1,label="PREDICTION\nCDS\nPRODIGAL",size=(100,100))
		button_prodigal.SetBackgroundColour("#CECECE")
		sizer.Add(button_prodigal,(1,1),(1,1),wx.EXPAND)

		button_protein=wx.Button(self,-1,label="PREDICTION\nPYL\nPROTEINS",size=(100,100))
		button_protein.SetBackgroundColour("#CECECE") 
		sizer.Add(button_protein,(1,2),(1,1),wx.EXPAND)

		button_blast=wx.Button(self,-1,label="BLAST",size=(100,100))
		button_blast.SetBackgroundColour("#CECECE") 
		sizer.Add(button_blast,(1,3),(1,1),wx.EXPAND) 

		sizer.Add(panel,(1,4),(1,1),wx.EXPAND)	

		button_pipeline=wx.Button(self,-1,label="COMPLETE PIPELINE")
		button_pipeline.SetBackgroundColour("lightgrey")
		font=wx.Font(12,wx.NORMAL,wx.NORMAL,wx.BOLD)
		button_pipeline.SetFont(font)
		sizer.Add(button_pipeline,(2,1),(1,3),wx.EXPAND)

		sizer.Add(panel,(3,0),(1,5),wx.EXPAND)

		self.Bind(wx.EVT_BUTTON,self.ouvrir_interface_protein,button_protein) #clic sur button_protein => lance ouvrir_interface_protein 
		self.Bind(wx.EVT_BUTTON,self.ouvrir_interface_blast,button_blast) 
		self.Bind(wx.EVT_BUTTON,self.ouvrir_interface_pipeline,button_pipeline)
		self.Bind(wx.EVT_BUTTON,self.ouvrir_interface_prodigal,button_prodigal)
		
	        self.SetSizerAndFit(sizer)
		self.Show(True)

	#ouverture interface prediction des proteines potentielles 
	def ouvrir_interface_protein(self,event): 
		home_directory=os.environ['HOME']
		cmd=home_directory+"/pyl_protein_prediction/Scripts/Interface/prediction_pyl.py"
		os.system(cmd) 

	#ouverture interface blast proteique + filtre résultats blast  
	def ouvrir_interface_blast(self,event):
		home_directory=os.environ['HOME']
		cmd=home_directory+"/pyl_protein_prediction/Scripts/Interface/blast_pyl.py"
		os.system(cmd) 

	#ouverture interface lancement pipeline complet 
	def ouvrir_interface_pipeline(self,event): 
		home_directory=os.environ['HOME']
		cmd=home_directory+"/pyl_protein_prediction/Scripts/Interface/pipeline.py"
		os.system(cmd)
	
	#ouverture interface lancement prodigal (prediction des CDS a partir d'un genome)
	def ouvrir_interface_prodigal(self,event):
		home_directory=os.environ['HOME'] 
		cmd=home_directory+"/pyl_protein_prediction/Scripts/Interface/prodigal.py" 
		os.system(cmd)
	       

if __name__ == "__main__":
    app = wx.App()
    frame = simpleapp_wx(None,-1,'Pipeline:Menu')
    app.MainLoop()
