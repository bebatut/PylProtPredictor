#!/usr/bin/env python
# -*- coding: iso-8859-1 -*-
# -*- coding: utf-8 -*-
import sys
import os
import wx
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)),
os.pardir))
import misc_functions
import predict_pyl_proteins


class simpleapp_wx(wx.Frame):
    def __init__(self, parent, id, title):
        wx.Frame.__init__(self, parent, id, title, size=(500, 350))
        self.parent = parent
        self.initialize()

    def initialize(self):
        sizer = wx.GridBagSizer(5, 15)

        vertical_space_20 = wx.Panel(self, -1, size=(0, 20))
        vertical_space_5 = wx.Panel(self, -1, size=(0, 5))

        title_font = wx.Font(13, wx.NORMAL, wx.NORMAL, wx.BOLD)
        validation_font = wx.Font(12, wx.NORMAL, wx.NORMAL, wx.BOLD)

        sizer.Add(vertical_space_20, (0, 0), (1, 5), wx.EXPAND)

        title = wx.StaticText(self, -1, label="Predict potential PYL proteins")
        title.SetFont(title_font)
        sizer.Add(title, (1, 1), (1, 3), wx.EXPAND)

        line = wx.StaticLine(self, -1)
        sizer.Add(line, (2, 0), (1, 5), wx.EXPAND)

        sizer.Add(vertical_space_5, (3, 0), (1, 5), wx.EXPAND)

        # sizer.Add(vertical_space_5, (5, 0), (1, 5), wx.EXPAND)

        genome_label = wx.StaticText(self, -1, label="Complete genome (fasta)")
        sizer.Add(genome_label, (6, 1), (1, 1), wx.EXPAND)
        genome_filepicker = wx.FilePickerCtrl(self, -1)
        sizer.Add(genome_filepicker, (7, 1), (1, 1), wx.EXPAND)

        predicted_cds_label = wx.StaticText(self, -1,
            label="Predicted CDS (fasta)")
        sizer.Add(predicted_cds_label, (6, 3), (1, 1), wx.EXPAND)
        predicted_cds_filepicker = wx.FilePickerCtrl(self, -1)
        sizer.Add(predicted_cds_filepicker, (7, 3), (1, 1), wx.EXPAND)

        # sizer.Add(vertical_space_5, (8, 0), (1, 5), wx.EXPAND)

        output_label = wx.StaticText(self, -1, label="Output directory")
        sizer.Add(output_label, (9, 1), (1, 1), wx.EXPAND)
        output_dirpicker = wx.DirPickerCtrl(self, -1)
        sizer.Add(output_dirpicker, (10, 1), (1, 1), wx.EXPAND)

        # sizer.Add(vertical_space_5, (11, 0), (1, 5), wx.EXPAND)

        validation_button = wx.Button(self, -1,
            label="Launch potential PYL protein prediction")
        validation_button.SetFont(validation_font)
        validation_button.SetBackgroundColour("lightgrey")
        sizer.Add(validation_button, (12, 1), (1, 3), wx.EXPAND)

        # sizer.Add(vertical_space_20, (13, 0), (1, 5), wx.EXPAND)

        self.Bind(wx.EVT_BUTTON, lambda event: self.launch_prediction(event,
            genome_filepicker.GetPath(),
            predicted_cds_filepicker.GetPath(),
            output_dirpicker.GetPath(), validation_button)

        self.SetSizerAndFit(sizer)
        self.Show(True)

    def launch_prediction(self, event, genome_filepath, predicted_cds_filepath,
        output_dirpath):
        if genome_filepath == "":
            wx.MessageBox('Missing file with complete genome',
                'Warning', wx.OK | wx.ICON_EXCLAMATION)
            return

        if not misc_functions.isfasta(genome_filepath):
            wx.MessageBox('Wrong format for file with complete genome',
                'Warning', wx.OK | wx.ICON_EXCLAMATION)
            return

        if predicted_cds_filepath == "":
            wx.MessageBox('Missing file with predicted CDS',
                'Warning', wx.OK | wx.ICON_EXCLAMATION)
            return

        if not misc_functions.isfasta(predicted_cds_filepath):
            wx.MessageBox('Wrong format for file with predicted CDS',
                'Warning', wx.OK | wx.ICON_EXCLAMATION)
            return

        if output_dirpath == "":
            wx.MessageBox('Missing output directory',
                'Warning', wx.OK | wx.ICON_EXCLAMATION)
            return

        predict_pyl_proteins.predict_pyl_proteins(
        genome_filepath, predicted_cds_filepath, output_dirpath)


if __name__ == "__main__":
    app = wx.App()
    frame = simpleapp_wx(None, -1,
        'Pipeline:prediction potential pyl proteins')
    app.MainLoop()
