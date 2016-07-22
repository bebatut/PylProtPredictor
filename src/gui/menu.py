#!/usr/bin/env python
# -*- coding: iso-8859-1 -*-
# -*- coding: utf-8 -*-
import os
import wx


class simpleapp_wx(wx. Frame):
    def __init__(self, parent, id, title):
        wx.Frame.__init__(self, parent, id, title, size=(400, 170))
        self.parent = parent
        self.initialize()

    def initialize(self):
        sizer = wx.GridBagSizer(20, 10)
        panel = wx.Panel(self, -1)
        sizer.Add(panel, (1, 0), (1, 1), wx.EXPAND)

        prodigal_button = wx.Button(self, -1,
            label="PREDICT\nCDS\nWITH\nPRODIGAL", size=(100, 100))
        prodigal_button.SetBackgroundColour("#CECECE")
        sizer.Add(prodigal_button, (1, 1), (1, 1), wx.EXPAND)
        self.Bind(wx.EVT_BUTTON, self.open_prodigal_gui, prodigal_button)

        pyl_protein_prediction_button = wx.Button(self, -1,
            label="PREDICT\nPYL\nPROTEINS", size=(100, 100))
        pyl_protein_prediction_button.SetBackgroundColour("#CECECE")
        sizer.Add(pyl_protein_prediction_button, (1, 2), (1, 1), wx.EXPAND)
        self.Bind(wx.EVT_BUTTON, self.open_pyl_protein_prediction_gui,
            pyl_protein_prediction_button)

        blast_button = wx.Button(self, -1, label="RUN\nBLAST", size=(100, 100))
        blast_button.SetBackgroundColour("#CECECE")
        sizer.Add(blast_button, (1, 3), (1, 1), wx.EXPAND)
        self.Bind(wx.EVT_BUTTON, self.open_blast_gui, blast_button)

        # sizer.Add(panel,(1,4),(1,1),wx.EXPAND)

        pipeline_button = wx.Button(self, -1, label="RUN WHOLE PIPELINE")
        pipeline_button.SetBackgroundColour("lightgrey")
        font = wx.Font(12, wx.NORMAL, wx.NORMAL, wx.BOLD)
        pipeline_button.SetFont(font)
        sizer.Add(pipeline_button, (2, 1), (1, 3), wx.EXPAND)
        self.Bind(wx.EVT_BUTTON, self.open_pipeline_gui, pipeline_button)

        # sizer.Add(panel, (3,0), (1,5), wx.EXPAND)

        self.SetSizerAndFit(sizer)
        self.Show(True)

    def open_pyl_protein_prediction_gui(self, event):
        os.system("pythonw src/gui/pyl_prediction_gui.py")

    def open_blast_gui(self, event):
        os.system("pythonw src/gui/blast_gui.py")

    def open_pipeline_gui(self, event):
        os.system("pythonw src/gui/pipeline_gui.py")

    def open_prodigal_gui(self, event):
        os.system("pythonw src/gui/prodigal_gui.py")


if __name__ == "__main__":
    app = wx.App()
    frame = simpleapp_wx(None, -1, 'Detection of pyrrolysine proteins: Menu')
    app.MainLoop()
