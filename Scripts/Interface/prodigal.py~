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

	#design+evenements
	def initialize(self):
	        sizer = wx.GridBagSizer(5,15)

		espace_vertical_20=wx.Panel(self,-1,size=(0,20))
		espace_vertical_5=wx.Panel(self,-1,size=(0,5))

		titre=wx.Font(13,wx.NORMAL,wx.NORMAL,wx.BOLD)		
		sous_titre=wx.Font(10,wx.NORMAL,wx.NORMAL,wx.BOLD)
		font_valider=wx.Font(12,wx.NORMAL,wx.NORMAL,wx.BOLD)
		
		sizer.Add(espace_vertical_20,(0,0),(1,5),wx.EXPAND)
		
		label_titre_prodigal=wx.StaticText(self,-1,label="PRODIGAL")
		label_titre_prodigal.SetFont(titre)
		sizer.Add(label_titre_prodigal,(1,1),(1,3),wx.EXPAND)	
		
		line=wx.StaticLine(self,-1)
		sizer.Add(line,(2,0),(1,5),wx.EXPAND)

		sizer.Add(espace_vertical_5,(3,0),(1,5),wx.EXPAND)

		label_genome=wx.StaticText(self,-1,label="Complete genome (fasta)") 
		sizer.Add(label_genome,(4,1),(1,1),wx.EXPAND) 
		filepicker_genome=wx.FilePickerCtrl(self,-1) 
		sizer.Add(filepicker_genome,(5,1),(1,1),wx.EXPAND) 

		label_result=wx.StaticText(self,-1,label="Results directory") 
		sizer.Add(label_result,(4,3),(1,1),wx.EXPAND)
		dirpicker_result=wx.DirPickerCtrl(self,-1) 
		sizer.Add(dirpicker_result,(5,3),(1,1),wx.EXPAND) 

		sizer.Add(espace_vertical_20,(6,0),(1,5),wx.EXPAND) 
		
		button_ok=wx.Button(self,-1,label="OK") 
		button_ok.SetBackgroundColour("lightgrey") 
		button_ok.SetFont(font_valider) 
		sizer.Add(button_ok,(7,1),(1,3),wx.EXPAND)
		
		sizer.Add(espace_vertical_20,(8,0),(1,5),wx.EXPAND)	

		#EVENEMENTS#
		self.Bind(wx.EVT_BUTTON,lambda event: self.launch_prodigal(event,filepicker_genome.GetPath(),dirpicker_result.GetPath()),button_ok)

		
		self.SetSizerAndFit(sizer)
		self.Show(True)	

	#lancement prediction CDS par Prodigal
	def launch_prodigal(self,event,path_genome,path_result):
		home_directory=os.environ['HOME']		
		path_prodigal=home_directory+"/pyl_protein_prediction/Prodigal"
		os.chdir(path_prodigal) 
		split_path_genome=path_genome.split("/") 
		file_genome=split_path_genome[len(split_path_genome)-1]
		name_genome=file_genome.split(".") 
		name_genome=name_genome[0]
		file_result=path_result+"/"+name_genome+"_predicted_CDS.fasta"
		cmd="./prodigal.linux -i "+path_genome+" -d "+file_result
		os.system(cmd)
		self.Destroy() 
		
		
		
			

if __name__ == "__main__":
    app = wx.App()
    frame = simpleapp_wx(None,-1,'Pipeline:Progigal')
    app.MainLoop()
