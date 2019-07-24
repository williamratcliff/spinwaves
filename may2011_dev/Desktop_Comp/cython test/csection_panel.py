"""

Disclaimer
==========

This software was developed at the National Institute of Standards and Technology at the NIST Center for 
Neutron Research by employees of the Federal Government in the course of their official duties. Pursuant 
to title 17 section 105* of the United States Code this software is not subject to copyright protection 
and is in the public domain. The SPINAL software package is an experimental spinwave analysis system. NIST 
assumes no responsibility whatsoever for its use, and makes no guarantees, expressed or implied, about its 
quality, reliability, or any other characteristic. The use of certain trade names or commercial products 
does not imply any endorsement of a particular product, nor does it imply that the named product is 
necessarily the best product for the stated purpose. We would appreciate acknowledgment if the software is 
used.

*Subject matter of copyright: United States Government works

Copyright protection under this title is not available for any work of the United States Government, but the 
United States Government is not precluded from receiving and holding copyrights transferred to it by 
assignment, bequest, or otherwise

Author: wflynn
"""

import wx
#import wxaddons.sized_controls as sc
import  wx.lib.intctrl
import  wx.grid as  gridlib
import numpy as np
import sys,os
import spinwaves.cross_section.csection_calc as csection_calc
import wx.richtext
from sympy import pi
import spinwaves.cross_section.util.printing as printing
from multiprocessing import Process, Pipe
import copy


class CrossSectionPanel(wx.Panel):
#    def __init__(self, procManager, *args, **kwds):
    def __init__(self, procManager, *args, **kwds):        
        kwds["style"] = wx.TAB_TRAVERSAL
        wx.Panel.__init__(self, *args, **kwds)
        
        self.tau_list = []
        
        self.window = wx.BoxSizer(wx.VERTICAL)
        
        self.file_box = wx.StaticBox(self, -1, "Select Files")
        self.int_file_label = wx.StaticText(self, -1, " Interaction File: ")
        self.int_browse_btn = wx.Button(self, -1, "Browse")
        self.int_file_txtCtrl = wx.TextCtrl(self, -1, "")
        self.spin_file_label = wx.StaticText(self, -1, " Spin File: ")
        self.spin_browse_btn = wx.Button(self, -1, "Browse")
        self.spin_file_txtCtrl = wx.TextCtrl(self, -1, "")
        self.tau_file_label = wx.StaticText(self, -1, " Tau File ")
        self.tau_browse_btn = wx.Button(self, -1, "Browse")
        self.tau_file_txtCtrl = wx.TextCtrl(self, -1, "")
        self.output_file_label = wx.StaticText(self, -1, " Output File ")
        self.output_browse_btn = wx.Button(self, -1, "Browse")
        self.output_file_txtCtrl = wx.TextCtrl(self, -1, "")

        self.hkl_interval = wx.StaticBox(self, -1, "HKL Interval")
        self.hklmin_label = wx.StaticText(self, -1, "hkl Min:")
        self.hklmin_txtCtrl = wx.TextCtrl(self, -1, "")
        self.hklmin_pi = wx.StaticText(self, -1, "*pi")
        self.hklmax_label = wx.StaticText(self, -1, "hkl Max:")
        self.hklmax_txtCtrl = wx.TextCtrl(self, -1, "")
        self.hklmax_pi = wx.StaticText(self, -1, "*pi")
        self.hkl_steps_label = wx.StaticText(self, -1, "Steps:")
        self.hkl_steps_ctrl = wx.SpinCtrl(self, -1, "", min=0, max=1000)
        
        self.w_interval = wx.StaticBox(self, -1, "Omega Interval")
        self.wmin_label = wx.StaticText(self, -1, "w Min:")
        self.wmin_txtCtrl = wx.TextCtrl(self, -1, "")
        self.wmax_label = wx.StaticText(self, -1, "w Max:")
        self.wmax_txtCtrl = wx.TextCtrl(self, -1, "")
        self.w_steps_label = wx.StaticText(self, -1, "Steps:")
        self.w_steps_ctrl = wx.SpinCtrl(self, -1, "", min=0, max=1000)
        
        self.scan_direction = wx.StaticBox(self, -1, "Scan Direction")
        self.kx_label = wx.StaticText(self, -1, "kx:")
        self.kx_txtCtrl = wx.TextCtrl(self, -1, "")
        self.ky_label = wx.StaticText(self, -1, "ky:")
        self.ky_txtCtrl = wx.TextCtrl(self, -1, "")
        self.kz_label = wx.StaticText(self, -1, "kz:")
        self.kz_txtCtrl = wx.TextCtrl(self, -1, "")
        self.temp_label = wx.StaticText(self, -1, "Temperature")
        self.temp_ctrl = wx.TextCtrl(self, -1, "")
        
        self.plot_range = wx.StaticBox(self, -1, "Plot Range")
        self.zmin_label = wx.StaticText(self, -1, "Minimum")
        self.zmin_ctrl = wx.TextCtrl(self, -1, "")
        self.zmax_label = wx.StaticText(self, -1, "Maximum")
        self.zmax_ctrl = wx.TextCtrl(self, -1, "")
        self.color_bar_box = wx.CheckBox(self, -1, "Display Color Bar")

        self.add_info = wx.StaticBox(self, -1, "Additional Information")
        self.spherical_avg_box = wx.CheckBox(self, -1, "Spherical Averaging")
        self.replot_box = wx.CheckBox(self, -1, "Just Re-plot?")
        self.ok_btn = wx.Button(self, -1, "Ok")
        self.cancel_btn = wx.Button(self, -1, "Cancel")        
        
        self.__set_properties()
        self.__do_layout()
        
        self.Bind(wx.EVT_BUTTON, self.OnIntFileBrowse, self.int_browse_btn)
        self.Bind(wx.EVT_BUTTON, self.OnSpinFileBrowse, self.spin_browse_btn)
        self.Bind(wx.EVT_BUTTON, self.OnTauFileBrowse, self.tau_browse_btn)
        self.Bind(wx.EVT_BUTTON, self.OnOutFileBrowse, self.output_browse_btn)
        self.Bind(wx.EVT_BUTTON, self.OnOk, self.ok_btn)
        self.Bind(wx.EVT_BUTTON, self.OnCancel, self.cancel_btn)

        self.processManager = procManager

    def __set_properties(self):
        # begin wxGlade: SpinwavePanel.__set_properties
        self.int_file_txtCtrl.SetMinSize((230, 27))
        self.spin_file_txtCtrl.SetMinSize((230, 27))
        self.tau_file_txtCtrl.SetMinSize((230, 27))
        self.output_file_txtCtrl.SetMinSize((230, 27))
        # end wxGlade

        self.hklmin_txtCtrl.SetValue(str(0))
        self.hklmax_txtCtrl.SetValue(str(2))
        self.hkl_steps_ctrl.SetValue(100)
        self.wmin_txtCtrl.SetValue(str(0))
        self.wmax_txtCtrl.SetValue(str(5))
        self.w_steps_ctrl.SetValue(100)
        
        self.kx_txtCtrl.SetValue(str(1))
        self.ky_txtCtrl.SetValue(str(0))
        self.kz_txtCtrl.SetValue(str(0))
        self.temp_ctrl.SetValue(str(0.0001))
        
        self.zmin_ctrl.SetValue(str(0))
        self.zmax_ctrl.SetValue(str(25))
        self.color_bar_box.SetValue(True)
        self.spherical_avg_box.SetValue(False)
        self.replot_box.SetValue(False)

    def __do_layout(self):
        # begin wxGlade: SpinwavePanel.__do_layout
       
        filebox = wx.StaticBoxSizer(self.file_box, wx.VERTICAL)
        intfile = wx.BoxSizer(wx.HORIZONTAL)
        spinfile = wx.BoxSizer(wx.HORIZONTAL)
        taufile = wx.BoxSizer(wx.HORIZONTAL)
        outfile = wx.BoxSizer(wx.HORIZONTAL)
        
        hklbox = wx.StaticBoxSizer(self.hkl_interval, wx.HORIZONTAL)
        hklmin = wx.BoxSizer(wx.VERTICAL)
        hklmax = wx.BoxSizer(wx.VERTICAL)
        hklsteps = wx.BoxSizer(wx.VERTICAL)
        hklminlab = wx.BoxSizer(wx.HORIZONTAL)
        hklminctrl = wx.BoxSizer(wx.HORIZONTAL)
        hklmaxlab = wx.BoxSizer(wx.HORIZONTAL)
        hklmaxtctrl = wx.BoxSizer(wx.HORIZONTAL)
        hklstepslab = wx.BoxSizer(wx.HORIZONTAL)
        hklstepsctrl = wx.BoxSizer(wx.HORIZONTAL)
        
        wbox = wx.StaticBoxSizer(self.w_interval, wx.HORIZONTAL)
        wmin = wx.BoxSizer(wx.VERTICAL)
        wmax = wx.BoxSizer(wx.VERTICAL)
        wsteps = wx.BoxSizer(wx.VERTICAL)
        wminlab = wx.BoxSizer(wx.HORIZONTAL)
        wminctrl = wx.BoxSizer(wx.HORIZONTAL)
        wmaxlab = wx.BoxSizer(wx.HORIZONTAL)
        wmaxtctrl = wx.BoxSizer(wx.HORIZONTAL)
        wstepslab = wx.BoxSizer(wx.HORIZONTAL)
        wstepsctrl = wx.BoxSizer(wx.HORIZONTAL)     
    
        scan = wx.StaticBoxSizer(self.scan_direction, wx.HORIZONTAL)
        kx = wx.BoxSizer(wx.VERTICAL)
        ky = wx.BoxSizer(wx.VERTICAL)
        kz = wx.BoxSizer(wx.VERTICAL)
        kxlab = wx.BoxSizer(wx.HORIZONTAL)
        kxctrl = wx.BoxSizer(wx.HORIZONTAL)
        kylab = wx.BoxSizer(wx.HORIZONTAL)
        kyctrl = wx.BoxSizer(wx.HORIZONTAL)
        kzlab = wx.BoxSizer(wx.HORIZONTAL)
        kzctrl = wx.BoxSizer(wx.HORIZONTAL)
        temp = wx.BoxSizer(wx.VERTICAL)
        templab = wx.BoxSizer(wx.HORIZONTAL)
        tempctrl = wx.BoxSizer(wx.HORIZONTAL)
        
        zmin = wx.BoxSizer(wx.VERTICAL)
        zmax = wx.BoxSizer(wx.VERTICAL)
        zminlab = wx.BoxSizer(wx.HORIZONTAL)
        zminctrl = wx.BoxSizer(wx.HORIZONTAL)
        zmaxlab = wx.BoxSizer(wx.HORIZONTAL)
        zmaxctrl = wx.BoxSizer(wx.HORIZONTAL)
        
        addinfo = wx.StaticBoxSizer(self.add_info, wx.HORIZONTAL)
        plotrange = wx.StaticBoxSizer(self.plot_range, wx.HORIZONTAL)
        buttons = wx.BoxSizer(wx.HORIZONTAL) 
        
        intfile.Add(self.int_file_label, 0, wx.ALIGN_CENTER_VERTICAL, 0)
        intfile.Add(self.int_file_txtCtrl, 0, wx.ALIGN_CENTER_VERTICAL, 0)
        intfile.Add(self.int_browse_btn, 0, wx.ALIGN_CENTER_VERTICAL, 0)
        filebox.Add(intfile, 0, wx.EXPAND, 0)
        spinfile.Add(self.spin_file_label, 0, wx.ALIGN_CENTER_VERTICAL, 0)
        spinfile.Add(self.spin_file_txtCtrl, 0, wx.ALIGN_CENTER_VERTICAL, 0)
        spinfile.Add(self.spin_browse_btn, 0, wx.ALIGN_CENTER_VERTICAL, 0)
        filebox.Add(spinfile, 0, wx.EXPAND, 0)
        taufile.Add(self.tau_file_label, 0 , wx.ALIGN_CENTER_VERTICAL, 0)
        taufile.Add(self.tau_file_txtCtrl, 0 , wx.ALIGN_CENTER_VERTICAL, 0)
        taufile.Add(self.tau_browse_btn, 0 , wx.ALIGN_CENTER_VERTICAL, 0)
        filebox.Add(taufile, 0, wx.EXPAND, 0)
        outfile.Add(self.output_file_label, 0 , wx.ALIGN_CENTER_VERTICAL, 0)
        outfile.Add(self.output_file_txtCtrl, 0 , wx.ALIGN_CENTER_VERTICAL, 0)
        outfile.Add(self.output_browse_btn, 0 , wx.ALIGN_CENTER_VERTICAL, 0)
        filebox.Add(outfile, 0, wx.EXPAND, 0)
        self.window.Add(filebox, 0, wx.EXPAND, 0)

        hklbox.Add((15, 15), 0, 0, 0)
        hklminlab.Add(self.hklmin_label, 0, wx.ALIGN_BOTTOM|wx.ALIGN_CENTER_HORIZONTAL, 0)
        hklmin.Add(hklminlab, 1, wx.EXPAND, 0)
        hklminctrl.Add(self.hklmin_txtCtrl, 0, wx.ALIGN_CENTER_HORIZONTAL|wx.ALIGN_CENTER_VERTICAL, 0)
        hklminctrl.Add(self.hklmin_pi, 0, wx.ALIGN_CENTER_VERTICAL, 0)
        hklmin.Add(hklminctrl, 1, wx.EXPAND, 0)
        hklbox.Add(hklmin, 1, wx.EXPAND, 0)
        hklmaxlab.Add(self.hklmax_label, 0, wx.ALIGN_BOTTOM|wx.ALIGN_CENTER_HORIZONTAL, 0)
        hklmax.Add(hklmaxlab, 1, wx.EXPAND, 0)
        hklmaxtctrl.Add(self.hklmax_txtCtrl, 0, wx.ALIGN_CENTER_HORIZONTAL|wx.ALIGN_CENTER_VERTICAL, 0)
        hklmaxtctrl.Add(self.hklmax_pi, 0, wx.ALIGN_CENTER_VERTICAL, 0)
        hklmax.Add(hklmaxtctrl, 1, wx.EXPAND, 0)
        hklbox.Add(hklmax, 1, wx.EXPAND, 0)
        hklstepslab.Add(self.hkl_steps_label, 0, wx.ALIGN_BOTTOM|wx.ALIGN_CENTER_HORIZONTAL, 0)
        hklsteps.Add(hklstepslab, 1, wx.EXPAND, 0)
        hklstepsctrl.Add(self.hkl_steps_ctrl, 0, wx.ALIGN_CENTER_HORIZONTAL|wx.ALIGN_CENTER_VERTICAL, 0)
        hklsteps.Add(hklstepsctrl, 1, wx.EXPAND, 0)
        hklbox.Add(hklsteps, 1, wx.EXPAND, 0)
        self.window.Add(hklbox, 0, wx.EXPAND, 0)

        wbox.Add((15, 15), 0, 0, 0)
        wminlab.Add(self.wmin_label, 0, wx.ALIGN_BOTTOM|wx.ALIGN_CENTER_HORIZONTAL, 0)
        wmin.Add(wminlab, 1, wx.EXPAND, 0)
        wminctrl.Add(self.wmin_txtCtrl, 0, wx.ALIGN_CENTER_HORIZONTAL|wx.ALIGN_CENTER_VERTICAL, 0)
        wmin.Add(wminctrl, 1, wx.EXPAND, 0)
        wbox.Add(wmin, 1, wx.EXPAND, 0)
        wmaxlab.Add(self.wmax_label, 0, wx.ALIGN_BOTTOM|wx.ALIGN_CENTER_HORIZONTAL, 0)
        wmax.Add(wmaxlab, 1, wx.EXPAND, 0)
        wmaxtctrl.Add(self.wmax_txtCtrl, 0, wx.ALIGN_CENTER_HORIZONTAL|wx.ALIGN_CENTER_VERTICAL, 0)
        wmax.Add(wmaxtctrl, 1, wx.EXPAND, 0)
        wbox.Add(wmax, 1, wx.EXPAND, 0)
        wstepslab.Add(self.w_steps_label, 0, wx.ALIGN_BOTTOM|wx.ALIGN_CENTER_HORIZONTAL, 0)
        wsteps.Add(wstepslab, 1, wx.EXPAND, 0)
        wstepsctrl.Add(self.w_steps_ctrl, 0, wx.ALIGN_CENTER_HORIZONTAL|wx.ALIGN_CENTER_VERTICAL, 0)
        wsteps.Add(wstepsctrl, 1, wx.EXPAND, 0)
        wbox.Add(wsteps, 1, wx.EXPAND, 0)
        self.window.Add(wbox, 0, wx.EXPAND, 0)

        scan.Add((15, 15), 0, 0, 0)
        kxlab.Add(self.kx_label, 0, wx.ALIGN_BOTTOM|wx.ALIGN_CENTER_HORIZONTAL, 0)
        kx.Add(kxlab, 1, wx.EXPAND, 0)
        kxctrl.Add(self.kx_txtCtrl, 0, wx.ALIGN_CENTER_HORIZONTAL|wx.ALIGN_CENTER_VERTICAL, 0)
        kx.Add(kxctrl, 1, wx.EXPAND, 0)        
        scan.Add(kx, 1, wx.EXPAND, 0)
        kylab.Add(self.ky_label, 0, wx.ALIGN_BOTTOM|wx.ALIGN_CENTER_HORIZONTAL, 0)
        ky.Add(kylab, 1, wx.EXPAND, 0)
        kyctrl.Add(self.ky_txtCtrl, 0, wx.ALIGN_CENTER_HORIZONTAL|wx.ALIGN_CENTER_VERTICAL, 0)
        ky.Add(kyctrl, 1, wx.EXPAND, 0)        
        scan.Add(ky, 1, wx.EXPAND, 0)
        kzlab.Add(self.kz_label, 0, wx.ALIGN_BOTTOM|wx.ALIGN_CENTER_HORIZONTAL, 0)
        kz.Add(kzlab, 1, wx.EXPAND, 0)
        kzctrl.Add(self.kz_txtCtrl, 0, wx.ALIGN_CENTER_HORIZONTAL|wx.ALIGN_CENTER_VERTICAL, 0)
        kz.Add(kzctrl, 1, wx.EXPAND, 0)        
        scan.Add(kz, 1, wx.EXPAND, 0)
        self.window.Add(scan, 0, wx.EXPAND, 0)
        
        plotrange.Add((15,15), 0, 0, 0)
        zminlab.Add(self.zmin_label, 0, wx.ALIGN_BOTTOM|wx.ALIGN_CENTER_HORIZONTAL, 0)
        zmin.Add(zminlab, 1, wx.EXPAND, 0)
        zminctrl.Add(self.zmin_ctrl, 0, wx.ALIGN_CENTER_HORIZONTAL|wx.ALIGN_CENTER_VERTICAL, 0)
        zmin.Add(zminctrl, 1, wx.EXPAND, 0)
        plotrange.Add(zmin, 1, wx.EXPAND, 0)
        zmaxlab.Add(self.zmax_label, 0, wx.ALIGN_BOTTOM|wx.ALIGN_CENTER_HORIZONTAL, 0)
        zmax.Add(zmaxlab, 1, wx.EXPAND, 0)
        zmaxctrl.Add(self.zmax_ctrl, 0, wx.ALIGN_CENTER_HORIZONTAL|wx.ALIGN_CENTER_VERTICAL, 0)
        zmax.Add(zmaxctrl, 1, wx.EXPAND, 0)
        plotrange.Add(zmax, 1, wx.EXPAND, 0)
        plotrange.Add(self.color_bar_box, 1, wx.ALIGN_BOTTOM|wx.ALIGN_CENTER_HORIZONTAL, 0)
        self.window.Add(plotrange, 0, wx.EXPAND, 0)

        addinfo.Add((15,15), 0, 0, 0)
        templab.Add(self.temp_label, 0, wx.ALIGN_BOTTOM|wx.ALIGN_CENTER_HORIZONTAL, 0)
        temp.Add(templab, 1, wx.EXPAND, 0)
        tempctrl.Add(self.temp_ctrl, 0, wx.ALIGN_CENTER_HORIZONTAL|wx.ALIGN_CENTER_VERTICAL, 0)
        temp.Add(tempctrl, 1, wx.EXPAND, 0)
        addinfo.Add(temp, 1, wx.EXPAND, 0)
        addinfo.Add(self.spherical_avg_box, 1, wx.ALIGN_BOTTOM|wx.ALIGN_CENTER_HORIZONTAL, 0)
        addinfo.Add(self.replot_box, 1, wx.ALIGN_BOTTOM|wx.ALIGN_CENTER_HORIZONTAL, 0)
        self.window.Add(addinfo, 0, wx.EXPAND, 0)

        buttons.Add(self.ok_btn, 0, wx.ALIGN_BOTTOM|wx.ALIGN_CENTER_HORIZONTAL, 0)
        buttons.Add(self.cancel_btn, 0, wx.ALIGN_BOTTOM|wx.ALIGN_CENTER_HORIZONTAL, 0)
        self.window.Add(buttons, 0, wx.ALIGN_CENTER_HORIZONTAL, 0)

        self.SetSizer(self.window)
        self.SetAutoLayout(True)
        self.window.Fit(self)
        self.SendSizeEvent()
        self.Refresh()
        # end wxGlade
        
    def OnIntFileBrowse(self, event):
        confBase = wx.ConfigBase.Create()
        confBase.SetStyle(wx.CONFIG_USE_LOCAL_FILE)
        defaultDir=confBase.Get().GetPath()

        wildcard="files (*.txt)|*.txt|All files (*.*)|*.*"
        int_dlg = wx.FileDialog(self, message="Choose an Interaction file", 
                            defaultDir=defaultDir, defaultFile="", wildcard=wildcard,
                            style=wx.OPEN | wx.CHANGE_DIR)

        if int_dlg.ShowModal() == wx.ID_OK:
            self.int_file_txtCtrl.SetValue(int_dlg.GetPath())
        int_dlg.Destroy()
    
    def OnSpinFileBrowse(self, event):
        confBase = wx.ConfigBase.Create()
        confBase.SetStyle(wx.CONFIG_USE_LOCAL_FILE)
        defaultDir=confBase.Get().GetPath()

        wildcard="files (*.txt)|*.txt|All files (*.*)|*.*"
        spin_dlg = wx.FileDialog(self, message="Choose a Spin Configuration file",
                            defaultDir=defaultDir, defaultFile="", wildcard=wildcard,
                            style=wx.OPEN | wx.CHANGE_DIR)

        if spin_dlg.ShowModal() == wx.ID_OK:
            self.spin_file_txtCtrl.SetValue(spin_dlg.GetPath())
        spin_dlg.Destroy()
    
    def OnTauFileBrowse(self, event):
        confBase = wx.ConfigBase.Create()
        confBase.SetStyle(wx.CONFIG_USE_LOCAL_FILE)
        defaultDir=confBase.Get().GetPath()

        wildcard="files (*.txt)|*.txt|All files (*.*)|*.*"
        tau_dlg = wx.FileDialog(self, message="Choose a file of tau vectors",
                            defaultDir=defaultDir, defaultFile="", wildcard=wildcard,
                            style=wx.OPEN | wx.CHANGE_DIR)

        if tau_dlg.ShowModal() == wx.ID_OK:
            self.tau_file_txtCtrl.SetValue(tau_dlg.GetPath())
        tau_dlg.Destroy()
    
    def OnOutFileBrowse(self, event):
        confBase = wx.ConfigBase.Create()
        confBase.SetStyle(wx.CONFIG_USE_LOCAL_FILE)
        defaultDir=confBase.Get().GetPath()

        wildcard="files (*.npz)|*.npz|All files (*.*)|*.*"
        out_dlg = wx.FileDialog(self, message="Choose a destination file for output",
                            defaultDir=defaultDir, defaultFile="", wildcard=wildcard,
                            style=wx.OPEN | wx.CHANGE_DIR)

        if out_dlg.ShowModal() == wx.ID_OK:
            self.output_file_txtCtrl.SetValue(out_dlg.GetPath())
        out_dlg.Destroy()
    
    def OnOk(self, event):
        failed, hkl_interval, w_interval, direction, temp, sphavg_bool, plotstats = self.Validate()
        
        if not failed:          
            int_file = self.int_file_txtCtrl.GetValue()
            spin_file = self.spin_file_txtCtrl.GetValue()
            output_file = self.output_file_txtCtrl.GetValue()
            
            tau_list = np.array(self.tau_list)
            
            self.processManager.startCrossSection(int_file, spin_file, hkl_interval, w_interval, tau_list, 
                                    direction, temp, output_file, sphavg_bool, plotstats)
#            csection_calc.cs_driver(int_file, spin_file, hkl_interval, w_interval, tau_list, 
#                                    direction, temp, output_file, sphavg_bool, plotstats)
#            self.processManager.startAnalyticDispersion(int_file, spin_file)
#            self.processManager.startNumericDispersion(int_file, spin_file, data, kMin*pi, kMax*p

    def OnCancel(self, event):
        self.GetParent().Close()
    
    def Validate(self):
        """Checks that all values are the right type. Any field that is not of the right
        type will be turned pink.
        
        Returns failed, data, kMin, kMax
        failed is True if validation fails and false otherwise."""
        
        hklmin = self.hklmin_txtCtrl.GetValue()
        hklmax = self.hklmax_txtCtrl.GetValue()
        hklsteps = self.hkl_steps_ctrl.GetValue()
        
        wmin = self.wmin_txtCtrl.GetValue()
        wmax = self.wmax_txtCtrl.GetValue()
        wsteps = self.w_steps_ctrl.GetValue()
        
        kx = self.kx_txtCtrl.GetValue()
        ky = self.ky_txtCtrl.GetValue()
        kz = self.kz_txtCtrl.GetValue()
        
        zmin = self.zmin_ctrl.GetValue()
        zmax = self.zmax_ctrl.GetValue()
        colorbar_bool = self.color_bar_box.GetValue()
        
        temp = self.temp_ctrl.GetValue()
        sphavg_bool = self.spherical_avg_box.GetValue()
        replot_bool = self.replot_box.GetValue()
         
        bgColor = "pink"
        failed = False
        
        #Validate hkl values
        num_hklmin = None
        num_hklmax = None
        try:
            num_hklmin = float(hklmin)*np.pi
            self.hklmin_txtCtrl.SetBackgroundColour("white")
        except:
            self.hklmin_txtCtrl.SetBackgroundColour(bgColor)
            failed = True
        try:
            num_hklmax = float(hklmax)*np.pi
            self.hklmax_txtCtrl.SetBackgroundColour("white")
        except:
            self.hklmax_txtCtrl.SetBackgroundColour(bgColor)
            failed = True      
        
        #Validate w values
        num_wmin = None
        num_wmax = None
        try:
            num_wmin = float(wmin)
            self.wmin_txtCtrl.SetBackgroundColour("white")
        except:
            self.wmin_txtCtrl.SetBackgroundColour(bgColor)
            failed = True
        try:
            num_wmax = float(wmax)
            self.wmax_txtCtrl.SetBackgroundColour("white")
        except:
            self.wmax_txtCtrl.SetBackgroundColour(bgColor)
            failed = True      
            
        #Validate kx,ky,kz,temp,zmin,zmax values
        num_kx = None
        num_ky = None
        num_kz = None
        num_temp = None
        num_zmin = None
        num_zmax = None
        try:
            num_kx = float(kx)
            self.kx_txtCtrl.SetBackgroundColour("white")
        except:
            self.kx_txtCtrl.SetBackgroundColour(bgColor)
            failed = True
        try:
            num_ky = float(ky)
            self.ky_txtCtrl.SetBackgroundColour("white")
        except:
            self.ky_txtCtrl.SetBackgroundColour(bgColor)
            failed = True      
        try:
            num_kz = float(kz)
            self.kz_txtCtrl.SetBackgroundColour("white")
        except:
            self.kz_txtCtrl.SetBackgroundColour(bgColor)
            failed = True
        try:
            num_temp = float(temp)
            self.temp_ctrl.SetBackgroundColour("white")
        except:
            self.temp_ctrl.SetBackgroundColour(bgColor)
            failed = True
        try:
            num_zmin = float(zmin)
            self.zmin_ctrl.SetBackgroundColour("white")
        except:
            self.zmin_ctrl.SetBackgroundColour(bgColor)
            failed = True
        try:
            num_zmax = float(zmax)
            self.zmax_ctrl.SetBackgroundColour("white")
        except:
            self.zmax_ctrl.SetBackgroundColour(bgColor)
            failed = True
            
        #Validate File Fields
        int_str = self.int_file_txtCtrl.GetValue()
        spin_str = self.spin_file_txtCtrl.GetValue()
        tau_str = self.tau_file_txtCtrl.GetValue()
        out_str = self.output_file_txtCtrl.GetValue()
        if int_str:
            self.int_file_txtCtrl.SetBackgroundColour("white")
        else: 
            self.int_file_txtCtrl.SetBackgroundColour(bgColor)
            failed = True
        if spin_str:
            self.spin_file_txtCtrl.SetBackgroundColour("white")
        else: 
            self.spin_file_txtCtrl.SetBackgroundColour(bgColor)
            failed = True
        if tau_str:
            self.tau_file_txtCtrl.SetBackgroundColour("white")
        else: 
            self.tau_file_txtCtrl.SetBackgroundColour(bgColor)
            failed = True
        if out_str:
            self.output_file_txtCtrl.SetBackgroundColour("white")
        else: 
            self.output_file_txtCtrl.SetBackgroundColour(bgColor)
            failed = True
            
        direction = {}
        direction['kx'] = num_kx
        direction['ky'] = num_ky
        direction['kz'] = num_kz
        hkl_interval = [num_hklmin, num_hklmax, int(self.hkl_steps_ctrl.GetValue())]
        w_interval = [num_wmin, num_wmax, int(self.w_steps_ctrl.GetValue())]
        
        tau_text = ''
        try:
            tau_file = open(tau_str,'r')
            tau_text = tau_file.read()
            self.tau_file_txtCtrl.SetBackgroundColour("white")
            print 'success', tau_text
        except:
            self.tau_file_txtCtrl.SetBackgroundColour(bgColor)
            failed = True
            print 'failed', tau_text

        items = tau_text.split()
        if len(items)%3 and not len(items):
            failed = True

        i = 0
        while not failed and i <= len(items)-3:
            tau1, tau2, tau3 = None, None, None
            try:
                tau1 = float(items[i])
                tau2 = float(items[i+1])
                tau3 = float(items[i+2])
                self.tau_file_txtCtrl.SetBackgroundColour("white")
            except:
                self.tau_file_txtCtrl.SetBackgroundColour(bgColor)
                failed = True
            self.tau_list.append([tau1,tau2,tau3])
            i+=3
        
        self.Refresh()
#        self.window.Show(True,True)
        
        plotstats = [num_zmin, num_zmax, colorbar_bool, replot_bool]
        
        return failed, hkl_interval, w_interval, direction, num_temp, sphavg_bool, plotstats


class MyApp(wx.App):
    def __init__(self, redirect=False, filename=None, useBestVisual=False, clearSigInt=True):
        wx.App.__init__(self,redirect,filename,clearSigInt)


    def OnInit(self):
        return True

def OnClose(event):
#    os.removedirs('spinwaves_temp')
    try:
        os.remove('csection_calc_data.npz')
    except:
        pass
    wx.Exit()

if __name__=='__main__':
    from spinwaves.utilities.Processes import ProcessManager
    app=MyApp()
    frame1 = wx.Frame(None, -1, "Spinwaves")#, size = (600,600))
    dlg=CrossSectionPanel(None, parent=frame1,id=-1)
    frame1.Bind(wx.EVT_CLOSE, OnClose)
    frame1.Show()

    app.MainLoop()
