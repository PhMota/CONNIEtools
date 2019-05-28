# coding: utf-8

import numpy as np

from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.backends.backend_gtk3cairo import (FigureCanvasGTK3Cairo as FigureCanvas)
from matplotlib.backends.backend_gtk3 import (NavigationToolbar2GTK3 as NavigationToolbar)
from matplotlib.figure import Figure
import matplotlib.ticker as ticker

import os
import time
import traceback

import gi
gi.require_version('Gtk', '3.0')
from gi.repository import Gtk, Gdk, GdkPixbuf, GObject, GLib
import threading

from ConnieDataPath import ConnieDataPath as ConniePaths
import ConnieImage
import Statistics

class ImageViewer(Gtk.Window):

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
        self.all_runIDs = ConniePaths.runID()
        self.update_runID = True
        self.runIDs = None
        
        self.add( self.build_window() )

        self.paths = []
        self.run = -1
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
        if self.run is not None:
            self.runEntry.get_child().set_text(str(self.run))

        self.connect( 'destroy', Gtk.main_quit )
        self.maximize()
        self.show_all()
        
        #if 'plot' in kwargs:
            #self.on_plotButton_click(None)
    
    def build_window( self ):
        body = Gtk.HBox()
        body.pack_start( self.build_optionsPanel(), False, False, 1 )
        body.pack_start( self.build_imagePanel(), True, True, 1 )
        body.pack_start( self.build_statisticsPanel(), False, False, 1 )
        return body

    def set_ohdu(self, ohdu ):
        self.ohduEntry.set_active( self.ohdus.index(ohdu) )
    
    def get_ohdu(self):
        try:
            return int( self.ohdus[ self.ohduEntry.get_active()] )
        except:
            return None
        
    def reset_axes(self):
        self.fig.clf()
        self.main_ax.cla()
        self.x_hist.cla()
        self.y_hist.cla()
        self.zoom_ax.cla()
        
    def build_imagePanel(self):
        self.fig = Figure(dpi=100, tight_layout=True)
        canvas = FigureCanvas(self.fig)
        self.fig.canvas.mpl_connect('draw_event', self.ondraw )
        toolbar = NavigationToolbar(canvas, self)
        children = toolbar.get_children()
        for i in range(len(children)-3):
            children[i].destroy()
        
        self.main_ax = None
        self.x_hist = None
        self.y_hist = None
        self.zoom_ax = None
        
        box = Gtk.VBox()
        box.pack_start(canvas, True, True, 0 )
        box.pack_end(toolbar, False, False, 0)
        return box

    def ondraw(self, event):
        pass
    
    def build_statisticsPanel(self):
        box = Gtk.VBox()
        self.zoom_fig = Figure()
        canvas = FigureCanvas( self.zoom_fig )
        b = Gtk.HBox()
        b.pack_start(canvas, True, True, 0)
        box.pack_start( b, True, True, 0)
        
        box.pack_start( Gtk.Label(label='Statistics', margin=10), False, False, 0)
        box.pack_start( Gtk.HSeparator(), False, False, 0)
        self.stats = Gtk.Label()
        self.stats.set_selectable(True)
        self.stats.set_use_markup(True)
        self.stats.set_justify(Gtk.Justification.RIGHT)
        self.stats.set_width_chars(40)
        
        scrolledStat = Gtk.ScrolledWindow()
        scrolledStat.set_policy(Gtk.PolicyType.NEVER, Gtk.PolicyType.AUTOMATIC)
        scrolledStat.add_with_viewport( self.stats )
        box.pack_end( scrolledStat, True, True, 0)
        
        return box
    
    def build_optionsPanel(self):
        optionsPanel = Gtk.VBox()
        
        optionsPanel.pack_start( Gtk.HBox(), False, False, 1 )
        optionsPanel.get_children()[-1].pack_start( Gtk.Label(label='run'), False, True, 1 )
        self.runEntry = Gtk.ComboBoxText( has_entry = True )
        optionsPanel.get_children()[-1].pack_end( self.runEntry, False, True, 1 )
        self.runEntry.get_child().set_width_chars(4)
        self.runList = []
        for run in ConniePaths.run()[::-1]:
            self.runEntry.append_text(str(run))
            self.runList.append(str(run))
        self.runEntry.connect('changed', self.on_run_changed)
        self.runEntry.get_child().connect( 'activate', self.on_plotButton_click )
        
        optionsPanel.pack_start( Gtk.HBox(), False, False, 1 )
        optionsPanel.get_children()[-1].pack_start( Gtk.Label(label='runID'), False, False, 1 )
        self.runIDEntry = Gtk.ComboBoxText( has_entry = True )
        self.runIDEntry.get_child().set_width_chars(5)
        self.runIDEntry.connect( 'changed', self.on_run_changed )
        self.runIDEntry.get_child().connect( 'activate', self.on_plotButton_click )
        optionsPanel.get_children()[-1].pack_end( self.runIDEntry, False, False, 1 )

        optionsPanel.pack_start( Gtk.HBox(), False, False, 1 )
        optionsPanel.get_children()[-1].pack_start( Gtk.Label(label='image'), False, False, 1 )
        self.imageTypeEntry = Gtk.ComboBoxText()
        self.imageTypeOptions = ['raw','osi','MB','mbs','scn']
        for key in self.imageTypeOptions:
            self.imageTypeEntry.append_text(key)
        optionsPanel.get_children()[-1].pack_end( self.imageTypeEntry, False, False, 1 )
        self.imageTypeEntry.connect( 'changed', self.on_plotButton_click )
        
        optionsPanel.pack_start( Gtk.HBox(), False, False, 1 )
        optionsPanel.get_children()[-1].pack_start( Gtk.Label(label='ohdu'), False, False, 1 )
        self.ohduEntry = Gtk.ComboBoxText()
        self.ohdus = [2,3,4,5,6,7,8,9,10,13,14]
        for ohdu in self.ohdus:
            self.ohduEntry.append_text(str(ohdu))
        optionsPanel.get_children()[-1].pack_end( self.ohduEntry, False, False, 1 )
        self.ohduEntry.connect( 'changed', self.on_plotButton_click )
        
        optionsPanel.pack_start( Gtk.HBox(), False, False, 1 )
        self.leftButton = Gtk.ToggleButton(label='left')
        optionsPanel.get_children()[-1].pack_start(self.leftButton, True, True, 0)
        self.leftButton.set_active(True)
        self.rightButton = Gtk.ToggleButton(label='right')
        optionsPanel.get_children()[-1].pack_end(self.rightButton, True, True, 0)
        self.leftButton.connect('toggled', self.on_side_toggled )
        self.rightButton.connect('toggled', self.on_side_toggled )
        
        self.plotOptions = ['E', 'x', 'y']
        self.optionEntry = {}
        for key in self.plotOptions:
            optionsPanel.pack_start( Gtk.HBox(), False, False, 1 )
            optionsPanel.get_children()[-1].pack_start( Gtk.Label(label=key), False, False, 1 )
            if key == 'x' or key =='y':
                for s in ['Max','Min']:
                    self.optionEntry['%s%s'%(key,s)] = Gtk.Entry()
                    self.optionEntry['%s%s'%(key,s)].set_text('auto')
                    self.optionEntry['%s%s'%(key,s)].set_width_chars(5)
                    optionsPanel.get_children()[-1].pack_end( self.optionEntry['%s%s'%(key,s)], False, False, 1 )
                    self.optionEntry['%s%s'%(key,s)].connect('activate', self.on_plotButton_click )
            else:
                self.optionEntry[key] = Gtk.Entry()
                optionsPanel.get_children()[-1].pack_end( self.optionEntry[key], False, False, 1 )
                self.optionEntry[key].connect('activate', self.on_plotButton_click )
        self.optionEntry['E'].set_text('200')
        self.optionEntry['E'].set_width_chars(7)

        self.plotButton = Gtk.Button(label='plot', relief=Gtk.ReliefStyle.NONE)
        #self.plotButton.set_relief()
        self.plotButton.connect( 'clicked', self.on_plotButton_click )
        optionsPanel.pack_end( self.plotButton, False, False, 1 )

        self.consoleButton = Gtk.ToggleButton(label='console')
        self.consoleButton.set_relief(Gtk.ReliefStyle.NONE)
        optionsPanel.pack_end( self.consoleButton, False, False, 1 )
        self.consolePopover = Gtk.Popover(modal=False)
        self.consolePopover.set_position(Gtk.PositionType.RIGHT)
        self.consolePopover.add( Gtk.VBox() )
        self.consolePopover.get_child().pack_start( Gtk.Label(label='Console', margin=10), False, False, 1)
        self.consolePopover.get_child().pack_start( Gtk.HSeparator(), False, False, 0 )
        self.Console = Gtk.Label(margin=10)
        self.consolePopover.get_child().pack_start( self.Console, True, True, 1)
        def toggled(widget):
            if widget.get_active():
                self.consolePopover.set_relative_to(widget)
                self.consolePopover.show_all()
            else:
                self.consolePopover.hide()
        self.consoleButton.connect('toggled', toggled)

        self.pathsLabel = Gtk.Label()
        self.pathsLabel.set_label( 'None' )
        self.pathsLabel.set_selectable(True)
        self.pathsLabel.set_property('margin', 10)
        self.pathsPopover = Gtk.Popover()
        self.pathsPopover.set_modal(False)
        self.pathsPopover.add( self.pathsLabel )
        pathsButton = Gtk.ToggleButton(label='paths')
        pathsButton.set_relief(Gtk.ReliefStyle.NONE)
        def toggled(widget):
            if widget.get_active():
                self.pathsPopover.set_relative_to(widget)
                self.pathsPopover.show_all()
            else:
                self.pathsPopover.hide()
        pathsButton.connect('toggled', toggled)
        optionsPanel.pack_end( pathsButton, False, False, 1 )

        return optionsPanel

    def refresh_runIDList(self, run=None, runID=None ):
        self.runIDEntry.remove_all()
        if run is None:
            self.runIDEntry.get_child().set_text('')
            return
        runIDs = ConniePaths.runID( run=run )[::-1]
        for runID_ in runIDs:
            self.runIDEntry.append_text(str(runID_))
        self.runIDs = runIDs
        if runID is None:
            runID = runIDs[0]
        print 'setting runID to', runID
        self.runIDEntry.get_child().set_text(str(runID))
        self.runIDEntry.set_active( self.runIDs.index(runID) )
        return

    def change_color(self, widget, color=None ):
        if color == 'red':
            color_ = Gdk.RGBA(red=1,green=0,blue=0)
        elif color is None:
            color_ = Gdk.RGBA(red=0,green=0,blue=0)
        else:
            print 'color not predicted', color
            return
        widget.override_color(0,color_)
        
    def on_run_changed( self, widget, *args ):
        if not self.update_runID:
            return
        self.update_runID = False
        if widget is self.runEntry:
            run = self.get_run()
            print 'on_run_changed', run
            if not str(run) in self.runList:
                print self.get_run(), 'is not present'
                self.change_color(self.runEntry.get_child(), 'red')
                self.refresh_runIDList()
            else:
                widget.set_active( self.runList.index(str(run)) )
                self.change_color( self.runEntry.get_child() )
                print 'run', run
                self.refresh_runIDList(run)
        elif widget is self.runIDEntry:
            runID = self.get_runID()
            print 'on_runID_changed', runID
            if not runID in self.all_runIDs:
                print self.get_runID(), 'is not present'
                self.change_color( self.runIDEntry.get_child(), 'red' )
                self.runEntry.get_child().set_text('')
            else:
                self.change_color(self.runIDEntry.get_child() )
                if runID not in self.runIDs:
                    print 'not present in self.runIDs'
                    run = ConniePaths.run(runID=runID)
                    self.runEntry.get_child().set_text(str(run))
                    self.runEntry.set_active( self.runList.index(str(run)) )
                    self.runEntry.get_child().override_color(0,Gdk.RGBA(red=0,green=0,blue=0))
                    print 'run', run
                    self.refresh_runIDList(run, runID)
                self.runIDEntry.set_active( self.runIDs.index(runID) )
        self.update_runID = True
        return
    
    def on_side_toggled(self, widget ):
        if widget.get_active():
            for child in widget.get_parent().get_children():
                if child is not widget: child.set_active(False)
            self.on_plotButton_click(widget)
        else:
            for child in widget.get_parent().get_children():
                if child is not widget: child.set_active(True)
        
    def on_runEntry_activate( self, entry ):
        self.run = int(entry.get_text())
        print 'activaterun', self.run
        self.runButtons[ self.run ].set_active(True)
    
    def on_runIDEntry_activate( self, entry ):
        print 'activaterunID', entry.get_text()
        self.runIDPopover.hide()
        self.on_plotButton_click(entry)
    
    def set_runID( self, runID ):
        #try:
            print 'set_runID', runID, type(runID)
            runID = int(runID)
            print 'set_runID datapath.run', ConniePaths.run(runID=runID)
            self.run = ConniePaths.run(runID=runID)
            print 'set_runID run', self.run
            #self.runButtons[self.run].set_active( True )
            self.runIDEntry.get_child().set_text(str(runID))
            #self.runIDEntry.set_active( self.runIDs.index(runID) )
            
            #self.runIDEntry.set_text( str(runID) )
            self.refresh_pathsLabel()
        #except:
            #return None

    def get_run( self ):
        try:
            return int(self.runEntry.get_child().get_text())
        except:
            return None

    def get_runID( self ):
        try:
            return int(self.runIDEntry.get_child().get_text())
        except:
            return None

    def get_sides( self ):
        return [ button.get_label() for button in sideButtons.get_children() if button.get_active() ]
    
    def set_imageType( self, imageType ):
        self.imageTypeEntry.set_active( self.imageTypeOptions.index(imageType) )
        return 
    
    def get_imageType( self ):
        return self.imageTypeOptions[ self.imageTypeEntry.get_active() ]
    
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

    def console_reset(self):
        self.Console.set_label('')
    
    def console_append(self, text ):
        label = self.Console.get_label()
        if label == '':
            self.Console.set_markup( '%s'%text )
        else:
            self.Console.set_markup( '%s\n%s'%(label, text) )
        while Gtk.events_pending(): Gtk.main_iteration()
        
    def on_plotButton_click( self, button ):
        self.console_reset()
        self.consoleButton.set_active(True)
        self.console_append('<span color="green">runID %s\nimageType %s\nohdu %s</span>'%(self.get_runID(), self.get_imageType(), self.get_ohdu() ))
        try:
            self.plotImage(button)
            self.console_append('<span color="green">plotting successful!</span>')
            self.consoleButton.set_active(False)
        except Exception as e:
            self.console_append('<span color="red">plotting failed :(</span>' )
            self.console_append('<span color="red">%s</span>'%str(e) )
            #self.console_append('<span color="red">%s</span>'%traceback.print_exc() )
            print traceback.print_exc()
        return

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
        
        if xMax <= xMin: 
            self.console_append('<span color="red">range_x [%s:%s]</span>'%(xMin,xMax))
            raise Exception('xMax &lt;= xMin')
        else: self.console_append('<span color="green">range_x [%s:%s]</span>'%(xMin,xMax))
        
        if yMax <= yMin: 
            self.console_append('<span color="red">range_y [%s:%s]</span>'%(yMin,yMax))
            raise Exception('yMax &lt;= yMin')
        else: self.console_append('<span color="green">range_y [%s:%s]</span>'%(yMin,yMax))
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
        ohdu = self.get_ohdu()
        
        if self.get_runID() is None:
            return False
    
        imageType = self.get_imageType()
        print self.get_runID(), self.get_ohdu(), imageType
        if imageType in ['osi','mbs']:
            imageType = 'raw'
        
        try:
            self.fullimage = ConnieImage.FullImage( runID = self.get_runID(), ohdu = self.get_ohdu(), imageType = imageType )
        except:
            raise Exception('%s not found'%imageType)
        self.console_append('<span color="green">image loaded</span>')

        if self.leftButton.get_active():
            self.image = self.fullimage.left()
            self.console_append('<span color="green">left side</span>')
        else:
            self.image = self.fullimage.right()
            self.console_append('<span color="green">right side</span>')

        if self.get_imageType() == 'osi':
            self.image = self.image.horizontalSubtraction()
            self.console_append('<span color="green">overscan subtracted</span>')
        
        if self.get_imageType() == 'mbs':
            self.image = self.image.horizontalSubtraction()
            self.console_append('<span color="green">overscan subtracted</span>')
            if self.rightButton.get_active():
                masterBias = ConnieImage.FullImage( runID = self.get_runID(), ohdu = self.get_ohdu(), imageType = 'MB' ).right()
            else:
                masterBias = ConnieImage.FullImage( runID = self.get_runID(), ohdu = self.get_ohdu(), imageType = 'MB' ).left()
            self.console_append('<span color="green">masterBias loaded</span>')
            self.image = ConnieImage.SideImage( self.image.image - masterBias.image )
            self.console_append('<span color="green">masterBias subtracted</span>')
        
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
        
        try:
            eRange = float(str(self.optionEntry['E'].get_text()))
            eMin = median - eRange
            eMax = median + eRange
            ebins = np.linspace( eMin, eMax, 100 )
        except:
            raise Exception('wrong energy range')

        self.main_ax.yaxis.set_major_formatter(ticker.FuncFormatter(lambda x, pos: '%s'%int(x+yMin)))
        self.main_ax.xaxis.set_major_formatter(ticker.FuncFormatter(lambda x, pos: '%s'%int(x+xMin)))
        section2 = section
        #if self.fast_render:
            #section2.image = section.image[::5,::5]
        section2.add_projection_to_axis( self.y_hist, axis=1, bins=ebins, align='vertical')
        self.fig.canvas.draw()
        self.console_append('<span color="green">X-projection done</span>')
        section2.add_projection_to_axis( self.x_hist, axis=0, bins=ebins, align='horizontal')
        self.fig.canvas.draw()
        self.console_append('<span color="green">Y-projection done</span>')
        section2.add_image_to_axis( self.main_ax, eMin, eMax )
        self.main_ax.yaxis.set_major_locator(ticker.AutoLocator())
        self.main_ax.xaxis.set_major_locator(ticker.AutoLocator())
        self.main_ax.yaxis.set_major_formatter(ticker.FuncFormatter(lambda x, pos: '%s'%int(x+yMin)))
        self.main_ax.xaxis.set_major_formatter(ticker.FuncFormatter(lambda x, pos: '%s'%int(x+xMin)))
        self.fig.canvas.draw()
        self.console_append('<span color="green">image done</span>')

        
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
            self.zoom_fig.canvas.draw()
            self.console_append('<span color="green">%s histogram done</span>'%sec.name )

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
        self.console_append('<span color="green">statistics done</span>')
        
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
 
