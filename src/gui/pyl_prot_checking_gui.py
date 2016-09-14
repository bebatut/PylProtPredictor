#!/usr/bin/env python
# -*- coding: iso-8859-1 -*-
# -*- coding: utf-8 -*-
import sys
import os
import re
import wx
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)),
os.pardir))
import check_pyl_protein
import misc_functions


class simpleapp_wx(wx.Frame):
    def __init__(self, parent, id, title):
        wx.Frame.__init__(self, parent, id, title)
        self.parent = parent
        self.initialize()

    def initialize(self):
        sizer = wx.GridBagSizer(5, 15)

        vertical_space_20 = wx.Panel(self, -1, size=(0, 20))
        vertical_space_5 = wx.Panel(self, -1, size=(0, 5))

        title_font = wx.Font(13, wx.NORMAL, wx.NORMAL, wx.BOLD)
        subtitle_font = wx.Font(10, wx.NORMAL, wx.NORMAL, wx.BOLD)
        validation_font = wx.Font(12, wx.NORMAL, wx.NORMAL, wx.BOLD)

        sizer.Add(vertical_space_20, (0, 0), (1, 5), wx.EXPAND)

        title_label = wx.StaticText(self, -1,
            label="Confirm potential PYL proteins with similarity search")
        title_label.SetFont(title_font)
        sizer.Add(title_label, (1, 1), (1, 5), wx.EXPAND)

        line = wx.StaticLine(self, -1)
        sizer.Add(line, (2, 0), (1, 5), wx.EXPAND)

        sizer.Add(vertical_space_5, (3, 0), (1, 5), wx.EXPAND)

        pot_pyl_prot_label = wx.StaticText(self, -1,
            label="Potential PYL protein to check (fasta file)")
        sizer.Add(pot_pyl_prot_label, (6, 1), (1, 1), wx.EXPAND)
        pot_pyl_prot_filepicker = wx.FilePickerCtrl(self, -1)
        sizer.Add(pot_pyl_prot_filepicker, (7, 1), (1, 1), wx.EXPAND)

        # sizer.Add(vertical_space_5, (9, 0), (1, 5), wx.EXPAND)

        ref_db_label = wx.StaticText(self, -1,
            label="Reference protein database (fasta file)")
        sizer.Add(ref_db_label, (8, 1), (1, 1), wx.EXPAND)
        ref_db_filepicker = wx.FilePickerCtrl(self, -1)
        sizer.Add(ref_db_filepicker, (9, 1), (1, 1), wx.EXPAND)

        evalue_label = wx.StaticText(self, -1,
            label="Maximum E-value (default, 10)")
        sizer.Add(evalue_label, (8, 3), (1, 1), wx.EXPAND)
        evalue = wx.TextCtrl(self, -1)
        sizer.Add(evalue, (9, 3), (1, 1), wx.EXPAND)

        # sizer.Add(vertical_space_5, (13, 0), (1, 5), wx.EXPAND)

        output_label = wx.StaticText(self, -1, label="Output directory")
        sizer.Add(output_label, (10, 1), (1, 1), wx.EXPAND)
        output_dirpicker = wx.DirPickerCtrl(self, -1)
        sizer.Add(output_dirpicker, (11, 1), (1, 1), wx.EXPAND)

        # sizer.Add(vertical_space_20, (17, 0), (1, 5), wx.EXPAND)

        launch_button = wx.Button(self, -1, label="Launch potential PYL protein checking")
        launch_button.SetBackgroundColour("lightgrey")
        launch_button.SetFont(validation_font)
        sizer.Add(launch_button, (12, 1), (1, 3), wx.EXPAND)

        # sizer.Add(vertical_space_20, (19, 0), (1, 5), wx.EXPAND)

        self.Bind(wx.EVT_BUTTON,
            lambda event: self.launch_pyl_protein_checking(event,
            pot_pyl_prot_filepicker.GetPath(),
            ref_db_filepicker.GetPath(),
            output_dirpicker.GetPath(),
            evalue.GetValue()), launch_button)

        self.SetSizerAndFit(sizer)
        self.Show(True)

    def launch_pyl_protein_checking(self, event, pot_pyl_prot_filepath, ref_db_filepath, output_dirpath, evalue):
        if pot_pyl_prot_filepath == "":
            wx.MessageBox('Missing file with the potential Pyl protein',
                'Warning', wx.OK | wx.ICON_EXCLAMATION)
            return

        if not misc_functions.isfasta(pot_pyl_prot_filepath):
            wx.MessageBox('Wrong format for file with the potential Pyl protein',
                'Warning', wx.OK | wx.ICON_EXCLAMATION)
            return

        if ref_db_filepath == "":
            wx.MessageBox('Missing file with reference database',
                'Warning', wx.OK | wx.ICON_EXCLAMATION)
            return

        if output_dirpath == "":
            wx.MessageBox('Missing output directory',
                'Warning', wx.OK | wx.ICON_EXCLAMATION)
            return

        if evalue == "":
            evalue = 10
        else:
            if not misc_functions.isfloat(evalue):
                wx.MessageBox('Evalue is not a float',
                    'Warning', wx.OK | wx.ICON_EXCLAMATION)
                return
            evalue = float(evalue)

        check_pyl_protein.check_potential_pyl_protein(pot_pyl_prot_filepath,
         ref_db_filepath, output_dirpath, evalue)

if __name__ == "__main__":
    app = wx.App()
    frame = simpleapp_wx(None, -1,
        'Check a potential Pyrrolysine protein')
    app.MainLoop()
