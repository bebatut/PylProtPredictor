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

		titre=wx.Font(13,wx.NORMAL,wx.NORMAL,wx.BOLD)		
		font_valider=wx.Font(12,wx.NORMAL,wx.NORMAL,wx.BOLD)

		sizer.Add(espace_vertical_20,(0,0),(1,5),wx.EXPAND)
		
		label_titre_prediction=wx.StaticText(self,-1,label="PREDICTION POTENTIAL PYL \nPROTEINS")
		label_titre_prediction.SetFont(titre)
		sizer.Add(label_titre_prediction,(1,1),(1,3),wx.EXPAND)	

		line=wx.StaticLine(self,-1)
		sizer.Add(line,(2,0),(1,5),wx.EXPAND)	

		sizer.Add(espace_vertical_5,(3,0),(1,5),wx.EXPAND)

		radio_button_scaf=wx.RadioBox(self,-1,choices=['Fully assembled genome','Scaffolds'])
		sizer.Add(radio_button_scaf,(4,1),(1,3),wx.EXPAND)

		sizer.Add(espace_vertical_5,(5,0),(1,5),wx.EXPAND)
	
		label_genome=wx.StaticText(self,-1,label="Complete genome (fasta)") 
		sizer.Add(label_genome,(6,1),(1,1),wx.EXPAND)

		filepicker_genome=wx.FilePickerCtrl(self,-1)
		sizer.Add(filepicker_genome,(7,1),(1,1),wx.EXPAND)

		label_cds=wx.StaticText(self,-1,label="Predicted CDS (fasta)")
		sizer.Add(label_cds,(6,3),(1,1),wx.EXPAND)

		filepicker_cds=wx.FilePickerCtrl(self,-1) 
		sizer.Add(filepicker_cds,(7,3),(1,1),wx.EXPAND)

		sizer.Add(espace_vertical_5,(8,0),(1,5),wx.EXPAND)

		label_result=wx.StaticText(self,-1,label="Results directory")
		sizer.Add(label_result,(9,1),(1,1),wx.EXPAND)

		dirpicker_result=wx.DirPickerCtrl(self,-1)
		sizer.Add(dirpicker_result,(10,1),(1,1),wx.EXPAND)

		sizer.Add(espace_vertical_20,(11,0),(1,5),wx.EXPAND)

		button=wx.Button(self,-1,label="OK") 
		button.SetFont(font_valider)
		button.SetBackgroundColour("lightgrey")
		sizer.Add(button,(12,1),(1,3),wx.EXPAND)

		sizer.Add(espace_vertical_20,(13,0),(1,5),wx.EXPAND)

		#EVENEMENTS#
		self.Bind(wx.EVT_BUTTON,lambda event: self.choice_script(event,filepicker_genome.GetPath(),filepicker_cds.GetPath(),dirpicker_result.GetPath(),radio_button_scaf.GetSelection()),button) 
		
		self.SetSizerAndFit(sizer)
		self.Show(True)

	def choice_script(self,event,path_genome,path_cds,path_result_dir,radio_box):
		os.system("clear") 		
		genome=True
		CDS=True
		results=True
		if path_genome=="": 
			print "Complete genome : select a file" 
			genome=False
		if path_cds=="": 
			print "Predicted CDS : select a file"
			CDS=False
		if path_result_dir=="": 
			print "Results directory : select a directory"
			results=False
		if genome and CDS and results:  
			if radio_box==0: #si "completely annoted genome" coché 
				self.launch_script(event,path_genome,path_cds,path_result_dir)
			elif radio_box==1: #si "scaffolds" coché 
				self.launch_script_scaffolds(event,path_genome,path_cds,path_result_dir)
	       
	def launch_script(self,event,path_genome,path_cds,path_result_dir): 
		home_directory=os.environ['HOME']
		script=home_directory+"/pyl_protein_prediction/Scripts/prediction_ORF_pyl_prodigal.py" 
		cmd="python "+script+" "+path_genome+" "+path_cds+" "+path_result_dir
		os.system(cmd) 
		
	def launch_script_scaffolds(self,event,path_genome,path_cds,path_result_dir): 
		home_directory=os.environ['HOME']
		script=home_directory+"/pyl_protein_prediction/Scripts/prediction_ORF_pyl_prodigal_scaffold.py" 
		cmd="python "+script+" "+path_genome+" "+path_cds+" "+path_result_dir
		os.system(cmd) 
		


if __name__ == "__main__":
    app = wx.App()
    frame = simpleapp_wx(None,-1,'Pipeline:prediction potential pyl proteins')
    app.MainLoop()
