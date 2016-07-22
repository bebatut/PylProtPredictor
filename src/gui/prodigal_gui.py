#!/usr/bin/env python
# -*- coding: iso-8859-1 -*-
# -*- coding: utf-8 -*-
import os
import wx
import sys
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)),
os.pardir))
import launch_prodigal
import misc_functions

class simpleapp_wx(wx.Frame):
    def __init__(self, parent, id, title):
        wx.Frame.__init__(self, parent, id, title, size=(500, 350))
        self.parent = parent
        self.initialize()

    def initialize(self):
        vertical_space_20 = wx.Panel(self, -1, size=(0, 20))
        vertical_space_5 = wx.Panel(self, -1, size=(0, 5))

        sizer = wx.GridBagSizer(5, 15)

        title_font = wx.Font(13, wx.NORMAL, wx.NORMAL, wx.BOLD)
        # subtitle_font = wx.Font(10, wx.NORMAL, wx.NORMAL, wx.BOLD)
        validation_font = wx.Font(12, wx.NORMAL, wx.NORMAL, wx.BOLD)

        # sizer.Add(vertical_space_20, (0, 0), (1, 5), wx.EXPAND)

        title = wx.StaticText(self, -1, label="Predict CDS with Prodigal")
        title.SetFont(title_font)
        sizer.Add(title, (1, 1), (1, 3), wx.EXPAND)

        line = wx.StaticLine(self, -1)
        sizer.Add(line, (2, 0), (1, 5), wx.EXPAND)

        sizer.Add(vertical_space_5, (3, 0), (1, 5), wx.EXPAND)

        genome_label = wx.StaticText(self, -1,
            label="Genome sequence (fasta)")
        sizer.Add(genome_label, (4, 1), (1, 1), wx.EXPAND)

        genome_filepicker = wx.FilePickerCtrl(self, -1)
        sizer.Add(genome_filepicker, (5, 1), (1, 1), wx.EXPAND)

        result_label = wx.StaticText(self, -1,
            label="Output directory")
        sizer.Add(result_label, (4, 3), (1, 1), wx.EXPAND)
        result_dirpicker = wx.DirPickerCtrl(self, -1)
        sizer.Add(result_dirpicker, (5, 3), (1, 1), wx.EXPAND)

        sizer.Add(vertical_space_20, (6, 0), (1, 5), wx.EXPAND)

        validation_button = wx.Button(self, -1, label="Launch Prodigal")
        validation_button.SetBackgroundColour("lightgrey")
        validation_button.SetFont(validation_font)
        sizer.Add(validation_button, (7, 1), (1, 3), wx.EXPAND)

        # sizer.Add(vertical_space_20, (8, 0), (1, 5), wx.EXPAND)

        self.Bind(wx.EVT_BUTTON, lambda event: self.launch_prodigal(event,
            genome_filepicker.GetPath(), result_dirpicker.GetPath()),
            validation_button)

        self.SetSizerAndFit(sizer)
        self.Show(True)

    def launch_prodigal(self, event, genome_filepath, result_dirpath):
        if genome_filepath == "":
            wx.MessageBox('Missing file with genome',
                'Warning', wx.OK | wx.ICON_EXCLAMATION)
            return

        if not misc_functions.isFasta(genome_filepath):
            wx.MessageBox('Wrong format for file with genome',
                'Warning', wx.OK | wx.ICON_EXCLAMATION)
            return

        if result_dirpath == "":
            wx.MessageBox('Missing output directory',
                'Warning', wx.OK | wx.ICON_EXCLAMATION)
            return

        launch_prodigal.launch_prodigal(genome_filepath, result_dirpath)
        self.Destroy()


if __name__ == "__main__":
    app = wx.App()
    frame = simpleapp_wx(None, -1,
        'Detection of pyrrolysine proteins: Progigal')
    app.MainLoop()
