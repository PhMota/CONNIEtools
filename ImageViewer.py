# coding: utf-8

import astropy
import astropy.io
import astropy.io.fits

import math
import numpy as np
from numpy.lib.recfunctions import append_fields, stack_arrays
import scipy
import scipy.stats
import scipy.signal
import scipy.special
from scipy.misc import factorial
import scipy.optimize

from collections import OrderedDict
import matplotlib
#matplotlib.use('Agg')
#matplotlib.use('qt4agg')
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from mpl_toolkits.axes_grid1 import make_axes_locatable, ImageGrid
#from matplotlib.backends.backend_gtk3agg import (FigureCanvasGTK3Agg as FigureCanvas)
from matplotlib.backends.backend_gtk3cairo import (FigureCanvasGTK3Cairo as FigureCanvas)
from matplotlib.backends.backend_gtk3 import (NavigationToolbar2GTK3 as NavigationToolbar)
from matplotlib.figure import Figure
import matplotlib.ticker as ticker

from functools import partial

import glob
import os
import shutil
import sys
import re
import copy
import root_numpy
import datetime
import time
import random
import traceback

import gi
gi.require_version('Gtk', '3.0')
from gi.repository import Gtk, Gdk, GdkPixbuf, GObject, GLib
import threading

from ConnieDataPath import ConnieDataPath as ConniePaths
import ConnieImage
import Statistics

class ImageViewer(Gtk.Window):
    #def __del__(self):
        #os.remove( self.tempPath )
        #Gtk.Window.__del__(self)

    def __init__(self, **kwargs ):
        self.id = '.' + os.environ['HOME'].split('/')[-1]
        self.id += str(int(time.time())%1000000)
        print 'id', self.id
        
        self.advanced = False
        if 'advanced' in kwargs:
            self.advanced = True
        
        self.tempPath = '/home/mota/public/gui/sessionimage%s.png'%self.id
        Gtk.Window.__init__( self, title="Image Viewer [session%s]"%self.id )
        self.maximize()

        self.set_border_width(3)
        self.image = None
        self.foot = None

        self.add( self.build_window() )

        self.paths = []
        self.run = -1
        self.runID_current = -1
        self.ohdu_current = -1
        self.imageType_current = -1
        self.no_file = 'file not found'
        
        self.fast_render = False
        if 'fast_render' in kwargs:
            self.fast_render = True
        if 'run' in kwargs: self.run = int(kwargs['run'])
        else: self.run = ConniePaths.run()[-1]
        if 'runID' in kwargs:
            self.set_runID( int(kwargs['runID']) )
            print 'runID', self.get_runID()
            print 'self.run', self.run
        if 'ohdu' in kwargs: self.set_ohdu( int(kwargs['ohdu']) )
        else: self.set_ohdu(2)
        self.runIDcolumn = None
        if 'imageType' in kwargs:
            self.set_imageType( kwargs['imageType'] )
        elif 'imagetype' in kwargs:
            self.set_imageType( kwargs['imagetype'] )
        else:
            self.set_imageType( 'raw' )
        
        print 'run', self.run
        if self.run is not None: self.runButtons[self.run].set_active(True)

        self.connect( 'destroy', Gtk.main_quit )
        self.maximize()
        self.show_all()
    
    def build_window(self):
        body = self.build_body()
        self.foot = self.build_foot()
        
        window = Gtk.VBox()
        window.pack_start( body, True, True, 3 )
        window.pack_start( self.foot, False, False, 3 )
        return window
    
    def build_foot(self):
        self.pathsLabel = Gtk.Label()
        self.pathsLabel.set_label( 'None' )
        self.pathsLabel.set_selectable(True)
        self.plotButton = Gtk.Button()
        self.plotButton.set_label('Plot')
        #self.plotButton.modify_bg( Gtk.StateType.NORMAL, Gdk.color_parse("blue") )
        self.plotButton.connect( 'clicked', self.on_plotButton_click )
        subbox = Gtk.VBox()
        subbox.pack_end( self.plotButton, False, False, 0 )
        foot = Gtk.HBox()
        foot.pack_end( subbox, False, False, 1 )
        foot.pack_start( self.pathsLabel, True, True, 1 )
        return foot

    def build_body( self ):
        runPanel = self.build_runPanel()
        runIDPanel = self.build_runIDPanel()
        mainPanel = self.build_mainPanel()
        
        body = Gtk.HBox()
        body.pack_start( runPanel, False, False, 1 )
        body.pack_start( runIDPanel, False, False, 1 )
        body.pack_end( mainPanel, True, True, 1 )
        return body

    def build_runPanel(self):
        runlabel = Gtk.Label()
        runlabel.set_label( 'runs' )

        runPanel = Gtk.VBox()
        runPanel.pack_start( runlabel, False, True, 1 )
        
        self.runEntry = Gtk.Entry()
        self.runEntry.set_width_chars(4)
        self.runEntry.connect('activate', self.on_runEntry_activate )
        runPanel.pack_start( self.runEntry, False, True, 1 )
        
        scrolledwindow = Gtk.ScrolledWindow()
        scrolledwindow.set_policy(Gtk.PolicyType.NEVER, Gtk.PolicyType.AUTOMATIC)
        runPanel.pack_start( scrolledwindow, True, True, 1 )

        self.runButtons = {}
        subbox = Gtk.VBox()
        for run in ConniePaths.run()[::-1]:
            self.runButtons[run] = Gtk.ToggleButton()
            self.runButtons[run].set_label( str(run) )
            self.runButtons[run].connect( 'toggled', self.on_runButton_toggle, run )
            subbox.pack_start( self.runButtons[run], False, False, 0 )
        scrolledwindow.add_with_viewport( subbox )

        return runPanel
        
    def build_runIDPanel(self):
        runIDcolumn = Gtk.VBox()

        self.runLabel = Gtk.Label()
        self.runLabel.set_label( 'runIDs' )
        runIDcolumn.pack_start( self.runLabel, False, False, 1 )

        self.runIDEntry = Gtk.Entry()
        self.runIDEntry.set_width_chars(5)
        self.runIDEntry.connect('activate', self.on_runIDEntry_activate )
        runIDcolumn.pack_start( self.runIDEntry, False, False, 1 )

        self.runIDScrolledWindow = Gtk.ScrolledWindow()
        self.runIDScrolledWindow.set_policy(Gtk.PolicyType.NEVER, Gtk.PolicyType.AUTOMATIC)
        runIDcolumn.pack_start( self.runIDScrolledWindow, True, True, 1 )
        
        return runIDcolumn

    def build_runIDButtons( self, run ):
        if len(self.runIDScrolledWindow.get_children()) > 0: self.runIDScrolledWindow.get_children()[0].destroy()
        runIDs = ConniePaths.runID(run=run)
        #print 'runIDs', runIDs
        subbox = Gtk.VBox()
        self.runIDButtons = {}
        for runID in runIDs[::-1]:
            self.runIDButtons[runID] = Gtk.ToggleButton()
            self.runIDButtons[runID].set_label( str(runID) )
            self.runIDButtons[runID].connect( 'clicked', self.on_runIDButton_toggle, runID )
            subbox.pack_start( self.runIDButtons[runID], False, False, 0 )
        self.runIDScrolledWindow.add_with_viewport(subbox)
        return
    
    def set_ohdu(self, ohdu ):
        self.optionEntry['ohdu'].set_text(str(ohdu))
    
    def get_ohdu(self):
        try:
            return int(self.optionEntry['ohdu'].get_text())
        except:
            return None
        
    def reset_axes(self):
        self.fig.clf()
        self.main_ax.cla()
        self.x_hist.cla()
        self.y_hist.cla()
        self.zoom_ax.cla()
        
    def build_figure(self):
        self.fig = Figure(dpi=100, tight_layout=True)

        self.canvas = FigureCanvas(self.fig)
        
        self.fig.canvas.mpl_connect('draw_event', self.ondraw )
        self.toolbar = NavigationToolbar(self.canvas, self)
        children = self.toolbar.get_children()
        for i in range(len(children)-3):
            children[i].destroy()
        
        self.main_ax = None
        self.x_hist = None
        self.y_hist = None
        self.zoom_ax = None
        #self.motion_notify_event = self.fig.canvas.callbacks['motion_notify_event']
        #self.fig.canvas.mpl_connect( 'motion_notify_event', self.on_mouse_move )
        return self.canvas, self.toolbar

    def ondraw(self, event):
        pass

    #def on_mouse_move( self, event ):
        #w, h = self.fig.canvas.get_width_height()
        #print float(event.x)/w, float(event.y)/h
        #if event.y/h > .9:
            #self.foot.set_visible(True)
        #else:
            #self.foot.set_visible(False)
            #self.foot.show()
        #self.toolbar.mouse_move(event)
        ##return
    
    def build_imagePanel(self):
        #self.imageCanvas = Gtk.Image()
        canvas, toolbar = self.build_figure()
        box = Gtk.VBox()
        box.pack_start(canvas, True, True, 0)
        box.pack_start(toolbar, False, False, 0)
        #scrolledwindowImage = Gtk.ScrolledWindow()
        #scrolledwindowImage.set_policy(Gtk.PolicyType.AUTOMATIC, Gtk.PolicyType.AUTOMATIC)
        #scrolledwindowImage.add_with_viewport( self.imageCanvas )
        #scrolledwindowImage.add_with_viewport( box )

        subbox = Gtk.VBox()
        #title = Gtk.Label()
        #title.set_text('Statistics')
        #subbox.pack_start(title, False, False, 3)
        
        self.zoom_fig = Figure()

        canvas = FigureCanvas( self.zoom_fig )
        b = Gtk.HBox()
        b.pack_start(canvas, True, True, 0)
        subbox.pack_start( b, True, True, 0)
        
        self.stats = Gtk.Label()
        self.stats.set_selectable(True)
        self.stats.set_use_markup(True)
        #self.stats.set_text()
        self.stats.set_justify(Gtk.Justification.RIGHT)
        self.stats.set_width_chars(40)
        
        scrolledStat = Gtk.ScrolledWindow()
        scrolledStat.set_policy(Gtk.PolicyType.NEVER, Gtk.PolicyType.AUTOMATIC)
        scrolledStat.add_with_viewport( self.stats )
        subbox.pack_end( scrolledStat, True, True, 0)
        
        imagePanel = Gtk.HBox()
        imagePanel.pack_start( box, True, True, 0 )
        imagePanel.pack_end( subbox, False, False, 0 )
        return imagePanel
    
    def build_mainPanel(self):
        optionsPanel = self.build_optionsPanel()
        imagePanel = self.build_imagePanel()
        
        mainPanel = Gtk.VBox()
        mainPanel.pack_start( optionsPanel, False, False, 1 )
        mainPanel.pack_start( imagePanel, True, True, 1 )
        return mainPanel
    
    def build_optionsPanel(self):
        optionsPanel = Gtk.VBox()
        firstLine = Gtk.HBox()

        self.plotOptions = ['ohdu','E', 'x', 'y']
        self.optionEntry = {}
        for key in self.plotOptions:
            optionLabel = Gtk.Label()
            optionLabel.set_label(key)
            firstLine.pack_start( optionLabel, False, False, 3 )
            if key == 'x' or key =='y':
                for s in ['Min','Max']:
                    self.optionEntry['%s%s'%(key,s)] = Gtk.Entry()
                    self.optionEntry['%s%s'%(key,s)].set_text('auto')
                    self.optionEntry['%s%s'%(key,s)].set_width_chars(5)
                    firstLine.pack_start( self.optionEntry['%s%s'%(key,s)], False, False, 0 )
                    #self.optionEntry['%s%s'%(key,s)].connect('activate', self.on_plotButton_click )
                    self.optionEntry['%s%s'%(key,s)].connect('activate', self.on_plotButton_click )
            else:
                self.optionEntry[key] = Gtk.Entry()
                firstLine.pack_start( self.optionEntry[key], False, False, 3 )
                #if key == 'E':
                    #self.optionEntry[key].connect('activate', self. )
                #else:
                    #self.optionEntry[key].connect('activate', self.on_plotButton_click )
                self.optionEntry[key].connect('activate', self.on_plotButton_click )
        self.optionEntry['E'].set_text('200')
        self.optionEntry['E'].set_width_chars(3)
        self.optionEntry['ohdu'].set_text('2')
        self.optionEntry['ohdu'].set_width_chars(2)
        
        self.Console = Gtk.Label()
        firstLine.pack_start(self.Console, True, True, 0)
        #self.sideButtons = Gtk.HBox()
        #firstLine.pack_end( self.sideButtons, False, False, 5 )

        #sideOptions = ['left', 'right']
        #for key in sideOptions:
            #button = Gtk.ToggleButton()
            #button.set_label( key )
            #self.sideButtons.pack_start( button, True, True, 0 )

        self.imageTypeButtons = Gtk.HBox()
        firstLine.pack_end( self.imageTypeButtons, False, False, 0 )
        #self.imageTypeOptions = ['raw','raw*','osi','mbs','scn','scn*']
        #imageTypeOptions = ['raw','osi','mbs','scn']
        imageTypeOptions = ['raw','osi','MB','mbs','scn']
        for key in imageTypeOptions:
            button = Gtk.ToggleButton()
            button.set_label( key )
            button.connect( 'toggled', self.on_imageTypeButton_toggle, key )
            self.imageTypeButtons.pack_start( button, True, True, 0 )

        #if self.advanced:
            #subtractOptions = ['R', '/R', '−R', '−MBR', '−MB', '−vOS', '−OS']
        #else:
            #subtractOptions = ['−vOS', '−OS']
        subtract = Gtk.HBox()
        subtractOptions = ['R']
        
        self.subtractButton = {}
        for key in subtractOptions:
            self.subtractButton[key] = Gtk.ToggleButton(key)
            subtract.pack_start( self.subtractButton[key], False, False, 0 )
            self.subtractButton[key].connect('toggled', self.on_plotButton_click )
        firstLine.pack_end( subtract, False, False, 10 )

        optionsPanel.pack_start(firstLine, False, True, 1)
        return optionsPanel

    def on_runEntry_activate( self, entry ):
        self.run = int(entry.get_text())
        print 'activaterun', self.run
        self.runButtons[ self.run ].set_active(True)
    
    def on_runIDEntry_activate( self, entry ):
        print 'activaterunID', entry.get_text()
        self.set_runID( int( entry.get_text() ) )
        self.on_plotButton_click(entry)
        #self.
    
    def set_runID( self, runID ):
        #try:
            print 'set_runID', runID, type(runID)
            runID = int(runID)
            print 'set_runID datapath.run', ConniePaths.run(runID=runID)
            self.run = ConniePaths.run(runID=runID)
            print 'set_runID run', self.run
            self.runButtons[self.run].set_active( True )
            if not runID in self.runIDButtons.keys(): self.build_runIDButtons( self.run )
            self.runIDButtons[ runID ].set_active( True )
            self.runIDEntry.set_text( str(runID) )
            self.refresh_pathsLabel()
        #except:
            #return None

    def get_runID( self ):
        try:
            return int(self.runIDEntry.get_text())
        except:
            return None

    def get_sides( self ):
        return [ button.get_label() for button in sideButtons.get_children() if button.get_active() ]
    
    def set_imageType( self, imageType ):
        for button in self.imageTypeButtons.get_children():
            if button.get_label() == imageType:
                button.set_active(True)

    def get_imageType( self ):
        for button in self.imageTypeButtons.get_children():
            if button.get_active() == True:
                return button.get_label()
        return None
    
    def refresh_pathsLabel( self ):
        runID = self.get_runID()
        imageType = self.get_imageType()
        self.plotButton.set_sensitive(False)
        print 'refresh', runID, imageType
        if runID is None or imageType is None: 
            print 'return'
            return
        paths = self.getPath()
        print paths
        self.pathsLabel.set_label( '\n'.join(paths) )
        if self.pathsLabel.get_text() != self.no_file: self.plotButton.set_sensitive(True)
    
    def get_paths(self):
        return self.pathsLabel.get_text().split('\n')
    
    def deactivateSiblings( self, button ):
        for sibling in button.get_parent().get_children():
            if sibling is not button:
                try:
                    sibling.set_active(False)
                except:
                    sibling.set_visited(False)
        return
            
    def on_imageTypeButton_toggle( self, button, key ):
        if button.get_active() == True:
            self.deactivateSiblings(button)
            self.refresh_pathsLabel()
            print 'toggle', key
            self.on_plotButton_click(None)
        
    def getPath( self ):
        runID = self.get_runID()
        imageType = self.get_imageType()
        print 'getPath', runID, imageType
        if runID is None or imageType is None: return [self.no_file]
        if imageType in ['raw','osi']:
            return ConniePaths.runIDPath( int(runID) )
        #elif imageType == 'osi':
            #return get_osi_path( int(runID) )
        elif imageType in ['MB','mbs']:
            return ConniePaths.masterBiasPath( int(runID) )
        elif imageType == 'scn':
            return ConniePaths.runIDPathProcessed( int(runID), image='scn' )
            #return get_scn_path( int(runID) )
        return [self.no_file]
    
    def on_runButton_toggle( self, button, run ):
        print 'toggle run', run, button.get_active()
        firstRunID = ConniePaths.runID( run=run )[0]
        if button.get_active() == True:
            self.deactivateSiblings( button )
            self.run = run
            self.build_runIDButtons( run )
            self.runEntry.set_text( str(self.run) )
        self.show_all()
        
    def on_runIDButton_toggle( self, button, runID ):
        print 'toggle runID', runID, button.get_active()
        if button.get_active() == True:
            self.deactivateSiblings(button)
            self.set_runID( runID )
            self.on_plotButton_click(None)
    
    def on_plotButton_click( self, button ):
        self.Console.set_markup('<span color="green">%s, runID%s[%s]</span>'%(self.get_imageType(), self.get_runID(), self.get_ohdu() ))
        while Gtk.events_pending(): Gtk.main_iteration()
        label = self.plotButton.get_label()
        #def running():
            #self.plotButton.set_label('running...')
            #self.plotButton.set_sensitive(False)
        #GLib.idle_add( running )
        self.plotButton.set_label('running...')
        self.plotButton.set_sensitive(False)
        try:
            self.plotImage(button)
            self.Console.set_markup('%s<span color="green">... Successful!</span>'%(self.Console.get_label()))
        except Exception as e:
            self.Console.set_markup('%s\n<span color="red">... Failed :([%s]</span>'%(self.Console.get_label(),str(e)) )
            print traceback.print_exc()
        
        self.plotButton.set_label(label)
        self.plotButton.set_sensitive(True)
        #def finished():
            #self.plotButton.set_label(label)
            #self.plotButton.set_sensitive(True)
        ##self.plotButton.show()
        #GLib.idle_add( finished )
        #while Gtk.events_pending(): Gtk.main_iteration()
    
    def parse_option( self, opt ):
        o = unicode(self.optionEntry[opt].get_text())
        return int(o) if o.isnumeric() else o
    
    def parse_lims( self ):
        xMin_str = self.optionEntry['xMin'].get_text()
        xMax_str = self.optionEntry['xMax'].get_text()
        yMin_str = self.optionEntry['yMin'].get_text()
        yMax_str = self.optionEntry['yMax'].get_text()
        
        xMin = None
        xMax = None
        yMin = None
        yMax = None

        if xMin_str == 'auto': xMin = 0
        if xMax_str == 'auto': xMax = self.image.width
        if xMin_str == 'os': xMin = -self.image.overscanWidth
        if xMax_str == 'os': xMax = -self.image.overscanWidth
        if yMin_str == 'auto': yMin = 0
        if yMax_str == 'auto': yMax = self.image.shape[0]
        if yMin_str == 'os': yMin = -self.image.overscanHeight
        if yMax_str == 'os': yMax = -self.image.overscanHeight
        
        if xMin is None: xMin = int(xMin_str)
        if yMin is None: yMin = int(yMin_str)
        if xMax is None: xMax = int(xMax_str)
        if yMax is None: yMax = int(yMax_str)
        
        if xMax < 0: xMax += self.image.width
        if yMax < 0: yMax += self.image.shape[0]
        
        if xMin < 0: xMin += xMax
        if yMin < 0: yMin += yMax
        
        if xMax_str.startswith('+'): xMax += xMin
        if yMax_str.startswith('+'): yMax += yMin
        
        print 'x,y', xMin,xMax,yMin,yMax
        return xMin,xMax,yMin,yMax
        
    def format_table( self, table, c ):
        lens = map( lambda col: max(map(len,col))+1, zip(*table) )
        result = ''
        for j in range(len(table[0])):
            result += ' '*(lens[j]-len(str(table[0][j])))+'<b>%s</b>'%(table[0][j])
        result +='\n'
        for i in range(1,len(table)):
            result += ' '*(lens[0]-len(str(table[i][0])))+'<span color="%s">%s</span>'%(c[i-1],table[i][0])
            for j in range(1,len(table[i])):
               result += ' '*(lens[j] - len(str(table[i][j])))+'<span>%s</span>'%table[i][j]
            result += '\n'
        #print 'result\n', result
        return result
            
        
    def plotImage( self, button ):
        stats = ''
        self.ohdu = int( self.optionEntry['ohdu'].get_text() )
        
        if self.get_runID() is None:
            return False
    
        imageType = self.get_imageType()
        print self.get_runID(), self.get_ohdu(), imageType
        if imageType in ['osi','mbs']:
            imageType = 'raw'
        
        self.fullimage = ConnieImage.FullImage( runID = self.get_runID(), ohdu = self.get_ohdu(), imageType = imageType )
        self.image = self.fullimage.left()

        if self.subtractButton['R'].get_active():
            self.image = self.fullimage.right()

        if self.get_imageType() == 'osi':
            self.image = self.image.horizontalSubtraction()
        
        if self.get_imageType() == 'mbs':
            self.image = self.image.horizontalSubtraction()
            if self.subtractButton['R'].get_active():
                masterBias = ConnieImage.FullImage( runID = self.get_runID(), ohdu = self.get_ohdu(), imageType = 'MB' ).right()
            else:
                masterBias = ConnieImage.FullImage( runID = self.get_runID(), ohdu = self.get_ohdu(), imageType = 'MB' ).left()
            self.image = ConnieImage.SideImage( self.image.image - masterBias.image )
        
        #if self.subtractButton['−R'].get_active():
            #B = np.mean(self.image.overscan().image)
            #print 'shapeB', B.shape
            #right = self.fullimage.right()
            #right.image += B - np.mean(right.overscan().image )
            #self.image = ConnieImage.SideImage( self.image.image - right.image )

        #if self.subtractButton['−OS'].get_active():
            #self.image = self.image.horizontalSubtraction()

        #if self.subtractButton['−vOS'].get_active():
            #self.image = self.image.verticalSubtraction()

        #if self.subtractButton['−MB'].get_active():
            #print 'masber Bias subtraction'
            #masterBias = ConnieImage.FullImage( runID = self.get_runID(), ohdu = self.get_ohdu(), imageType = 'MB' ).left()
            #cov = np.mean( np.median(masterBias.image,axis=0) * np.median(self.image.image, axis=0) )
            #a = cov/np.var( np.median(masterBias.image,axis=0))
            #print 'a =', a
            #self.image.image -= masterBias.image

        #if self.subtractButton['−MBR'].get_active():
            #masterBias = ConnieImage.FullImage( runID = self.get_runID(), ohdu = self.get_ohdu(), imageType = 'MB' ).right()
            #cov = np.mean( np.median(masterBias.image,axis=0) * np.median(self.image.image, axis=0) )
            #a = cov/np.var( np.median(masterBias.image,axis=0) )
            #print 'a =', a
            #self.image.image -= a*masterBias.image

        #if self.advanced:
            #if self.subtractButton['−OS'].get_active():
                #right = self.fullimage.right()
                #if self.subtractButton['−OS'].get_active():
                    #right = right.horizontalSubtraction()
                #if self.subtractButton['−vOS'].get_active():
                    #right = right.verticalSubtraction()
                #R = right.overscan().image - np.median(right.overscan().image, axis=1)[:,None]
                #L = self.image.overscan().image - np.median(self.image.overscan().image, axis=1)[:,None]
                #covLR = np.mean( L*R )
                #a = covLR/np.var(R)
                #b = Statistics.MAD(L)/Statistics.MAD(R)
                #print 'covLR', covLR, a, b
                #self.image = ConnieImage.SideImage( self.image.image - b*right.image )
            #elif self.subtractButton['−vOS'].get_active():
                #right = self.fullimage.right()
                #if self.subtractButton['−OS'].get_active():
                    #right = right.horizontalSubtraction()
                #if self.subtractButton['−vOS'].get_active():
                    #right = right.verticalSubtraction()
                #right = self.fullimage.right()
                #R = right.active().halfverticalOverScan().image - np.median(right.active().halfverticalOverScan().image, axis=1)[:,None]
                #L = self.image.active().halfverticalOverScan().image - np.median(self.image.active().halfverticalOverScan().image, axis=1)[:,None]
                #a = np.mean(L*R)/np.var(R)
                #print 'covLR_vos', np.mean(L*R), a
                #self.image = ConnieImage.SideImage( self.image.image - a*right.image )

        
        stats += 'shape (%s, %s)'%(self.image.shape[0],self.image.shape[1])
        stats += '\nrunID %s ohdu %s %s'%(self.get_runID(), self.get_ohdu(), self.get_imageType())

        if self.main_ax is not None:
            self.reset_axes()
        
        xMin,xMax,yMin,yMax = self.parse_lims()
        section = ConnieImage.ImageBase( self.image.image[yMin:yMax,xMin:xMax][::-1,:] )
        section.name = 'selected'

        self.main_ax = self.fig.add_subplot(111)
        divider = make_axes_locatable(self.main_ax)
        self.x_hist = divider.append_axes('top', size=1.2, pad=.1, sharex=self.main_ax)
        self.y_hist = divider.append_axes('right', size=1.2, pad=.1, sharey=self.main_ax)

        self.fig.canvas.draw()

        self.zoom_ax = self.zoom_fig.add_subplot(111)
        self.zoom_ax.set_yscale('log')
        self.zoom_ax.grid( True )
        self.zoom_ax.yaxis.set_label_position('right')
        self.zoom_ax.xaxis.set_label_position('top')
        self.zoom_ax.set_xlabel('E[adu]')
        self.zoom_ax.set_ylabel('counts')

        self.zoom_fig.canvas.draw()

        median = self.image.median()
        mad = self.image.MAD()
        mid = lambda _: .5*(_[1:]+_[:-1])
        diff = lambda _: _[1:]-_[:-1]
        
        eRange = float(str(self.optionEntry['E'].get_text()))
        eMin = median - eRange
        eMax = median + eRange
        ebins = np.linspace( eMin, eMax, 100 )

        self.main_ax.yaxis.set_major_formatter(ticker.FuncFormatter(lambda x, pos: '%s'%int(x+yMin)))
        self.main_ax.xaxis.set_major_formatter(ticker.FuncFormatter(lambda x, pos: '%s'%int(x+xMin)))
        section2 = section
        #if self.fast_render:
            #section2.image = section.image[::5,::5]
        section2.add_projection_to_axis( self.y_hist, axis=1, bins=ebins, align='vertical')
        self.fig.canvas.draw()
        section2.add_projection_to_axis( self.x_hist, axis=0, bins=ebins, align='horizontal')
        self.fig.canvas.draw()
        section2.add_image_to_axis( self.main_ax, eMin, eMax )
        self.main_ax.yaxis.set_major_locator(ticker.AutoLocator())
        self.main_ax.xaxis.set_major_locator(ticker.AutoLocator())
        self.main_ax.yaxis.set_major_formatter(ticker.FuncFormatter(lambda x, pos: '%s'%int(x+yMin)))
        self.main_ax.xaxis.set_major_formatter(ticker.FuncFormatter(lambda x, pos: '%s'%int(x+xMin)))
        self.fig.canvas.draw()

        
        table = [ ['section', 'shape', 'med', 'MAD'] ]
        cs = ['red', 'green', 'blue', 'yellow', 'magenta', 'black']
        median = {}
        mad = {}
        nohits = self.image.active().data().void( thr=15 )
        if nohits is not None:
            nohits.name = 'nohits'

        for c, sec in zip( cs, [ self.image.active().halfdata(), self.image.overscan(), self.image.active().halfverticalOverScan(), section, nohits ] ):
            if sec is None: continue
            print sec.name, 
            hist = self.zoom_ax.hist( sec.flatten(), bins=ebins, histtype='step', color=c, label='L' )[0]
            if len(sec.shape)>1:
                median[sec.name] = sec.mean_medians()
                mad[sec.name] = sec.mean_MADs()
                size = '% 4d,% 4d'%(sec.shape[0],sec.shape[1])
            else:
                median[sec.name] = sec.median()
                mad[sec.name] = sec.MAD()
                size = '%d'%sec.size
            table.append( [sec.name, size, '%.5g'%median[sec.name], '%.5g'%mad[sec.name]] )
        print
        stats += '\n' + self.format_table(table,cs)
        
        table = [ ['', 'os', 'vos'] ]
        if nohits is not None:
            Dmed_os = float(median['nohits'] - median['overscan'])
            Dmed_vos = float(median['nohits'] - median['hvertOS'])
            table.append( [ u'dE', '%.5g'%Dmed_os, '%.5g'%Dmed_vos ] )
            DVar_os = mad['nohits']**2 - mad['overscan']**2
            DVar_vos = mad['nohits']**2 - mad['hvertOS']**2
            table.append( [ u'dVar', '%.5g'%DVar_os, '%.5g'%DVar_vos ] )
            gain_os = DVar_os**2/Dmed_os**3
            gain_vos = DVar_vos**2/Dmed_vos**3
            if self.advanced:
                table.append( [ u'gain', '%.5g'%gain_os, '%.5g'%gain_vos ] )
                lamb_os = Dmed_os**2/DVar_os
                lamb_vos = Dmed_vos**2/DVar_vos
                table.append( [ 'lambda', '%.5g'%lamb_os, '%.5g'%lamb_vos ] )
                table.append( [ u'sigma', '%.5g'%(mad['overscan']/gain_os), '%.5g'%(mad['hvertOS']/gain_vos) ] )
        else:
            Dmed_os = float(median['hdata'] - median['overscan'])
            Dmed_vos = float(median['hdata'] - median['hvertOS'])
            table.append( [ u'dE', '%.5g'%Dmed_os, '%.5g'%Dmed_vos ] )
            DVar_os = mad['hdata']**2 - mad['overscan']**2
            DVar_vos = mad['hdata']**2 - mad['hvertOS']**2
            table.append( [ u'dVar', '%.5g'%DVar_os, '%.5g'%DVar_vos ] )
        stats += self.format_table(table, ['black']*10)
        
        self.stats.set_markup( stats + '\n' + self.stats.get_label()  )

        self.fig.canvas.draw()
        a, b = self.zoom_ax.get_ylim()
        self.zoom_ax.set_ylim((min(.5,a), b))
        self.zoom_fig.canvas.draw()
        
        print 'image ready'
        return True
 
