# coding: utf-8

import numpy as np

import matplotlib
import matplotlib.pyplot as plt
from matplotlib import patches, colors
from mpl_toolkits.axes_grid1 import make_axes_locatable, ImageGrid
from matplotlib.backends.backend_gtk3cairo import (FigureCanvasGTK3Cairo as FigureCanvas)
from matplotlib.backends.backend_gtk3 import (NavigationToolbar2GTK3 as NavigationToolbar)
from matplotlib.figure import Figure
import matplotlib.ticker as ticker

import glob
import os
import sys
import shutil
import time
import re

import gi
gi.require_version('Gtk', '3.0')
from gi.repository import Gtk, Gdk, GdkPixbuf, GObject, GLib
import threading

from ConnieDataPath import ConnieDataPath as ConniePaths
import ConnieImage
import Statistics as stats
import ImageViewer

class utils:
    @staticmethod
    def try_convertion( s, t ):
        result = None
        try: result = t(s)
        except: pass
        return result
    
def file_age(filepath):
    return time.time() - os.path.getmtime(filepath)

def file_age_str(filepath):
    secs = file_age(filepath)
    minutes = int(secs)/60
    hours = minutes/60
    return '%sh%sm'%(hours,minutes)

def open_imageViewer( runID, ohdu ):
    win = ImageViewer.ImageViewer( runID=runID, ohdu=ohdu, imageType='osi', plot=True )
    Gtk.main()

class Logger(object):
    def __init__( self, session, out, key ):
        self.terminal = out
        self.log_file = '/home/mota/public/gui/%s.%s'%(session,key)
        self.key = key
        self.session = session
        self.queue = []
        self.locked = False
    
    def write( self, message ):
        if self.locked: self.queue.append(message)
        self.locked = True
        log_message = '%s/%s[%s]: %s'%(self.session, self.key, time.strftime("%Y%b%d,%H:%M:%S"), message)
        self.terminal.write(message)
        self.terminal.flush()
        f = open( self.log_file, 'a+', 0 )
        f.write(log_message)
        f.close()
        self.locked = False
        if len(self.queue) > 0:
            self.write( '[queued]%s'%self.queue.pop(0) )
        return

class MonitorViewer(Gtk.Window):
    
    def __del__(self, widget=None, msg=None ):
        self.destroy()
        Gtk.main_quit(widget)
        print 'safely quitting...'
        self.quit = True
        #os.remove( self.tempPath )
        #self.remove_darkCurrentTable_lock()
        #self.remove_lock( self.darkCurrentRawTable )
        #self.remove_lock( self.readoutNoiseRawTable )
        return True
        
    def __init__(self):
        os.umask(0)

        self.id = '.' + os.environ['HOME'].split('/')[-1]
        self.id += str(int(time.time())%1000000)
        
        print 'id %s'%self.id
        sys.stdout = Logger( self.id, sys.stdout, 'out')
        sys.stderr = Logger( self.id, sys.stderr, 'err')
        
        self.log_file = '/home/mota/public/gui/%s.log'%self.id
        self.quantities = ['darkCurrent', 'readoutNoise','diffMADSqr']
        
        self.ohdus = range(2,16)
        self.global_runIDMax = None
        self.global_runIDMin = None
        self.runRanges = None
        
        self.firstLockMsg = {}
        self.firstAllDoneMsg = {}
        
        self.imagePaths = {}
        self.tablePaths = {}
        self.lockTimers = {}
        self.has_remove_lock_button = {}
        self.percentage_done = {}
        self.plotting = {}
        self.wait = {}
        self.thread = {}
        for quantity in self.quantities:
            self.imagePaths[quantity] = '/home/mota/public/gui/%sRawImage%s.png'%(quantity,self.id)
            self.tablePaths[quantity] = '/home/mota/public/gui/%sRawTable.csv'%(quantity)
            self.lockTimers[quantity] = None
            self.has_remove_lock_button[quantity] = False
            self.percentage_done[quantity] = None
            self.plotting[quantity] = False
            self.wait[quantity] = 0
            self.thread[quantity] = None
            self.firstAllDoneMsg[quantity] = True
            self.firstLockMsg[quantity] = True

        Gtk.Window.__init__( self, title="Monitor Viewer [session%s]"%self.id )
        self.connect( 'delete-event', self.__del__ )
        self.maximize()
        self.set_border_width(3)
        self.add( self.build_window() )
        print 'window built successfully'
        
        self.update = True
        print 'first plottig, number of threads %s'%threading.active_count()
        #for quantity in self.quantities:
            #self.plot_table( quantity )
        self.plot_table( self.currentPageLabel )

        #self.connect( 'motion-notify-event', self.on_mouse_move )
        self.show_all()
        self.quit = False

        for quantity in self.quantities:
            print 'request thread for %s %s'%(quantity, threading.active_count())
            self.start_thread( lambda quantity=quantity: self.computation_loop(quantity), quantity )

    def computation_loop( self, quantity ):
        print 'start %s loop in thread %s number of threads %s'%( quantity, threading.currentThread().getName(), threading.active_count() )
        while 1:
            if self.wait[quantity] > 0:
                #print 'sleeping', quantity, self.wait[quantity]
                for i in range(self.wait[quantity]):
                    time.sleep(1)
                    if self.should_quit():
                        exit(0)
                self.wait[quantity] = 0
            self.update_table( quantity )
            if self.should_quit():
                exit(0)
        return True
        
    def build_window(self):
        vbox = Gtk.VBox()
        header = self.build_header()
        body = self.build_body()
        vbox.pack_start(header, False, False, 1)
        vbox.pack_start(body, True, True, 1)
        return vbox
    
    def build_header(self):
        subbox = Gtk.HBox()
        label = Gtk.Label()
        label.set_label('Dark Current and Readout Noise')
        runIDLabel = Gtk.Label()
        runIDLabel.set_label('runIDs')
        self.runIDRangeEntry = Gtk.Entry()
        self.runIDRangeEntry.set_text('-100:auto')
        self.runIDRangeEntry.set_width_chars(len(' auto: auto'))
        self.runIDRangeEntry.connect( 'activate', self.on_runIDRangeEntry_activate )
        self.omitOutliersButton = Gtk.ToggleButton()
        self.omitOutliersButton.set_label('outliers')
        self.omitOutliersButton.set_active(True)
        self.omitOutliersButton.connect( 'clicked', self.on_omitOutliersButton_clicked )
        subbox.pack_start(label, True, True, 1 )
        subbox.pack_start(runIDLabel, False, False, 1 )
        subbox.pack_start(self.runIDRangeEntry, False, False, 1 )
        subbox.pack_start(self.omitOutliersButton, False, False, 1 )
        
        self.ohdusBox = Gtk.HBox()
        subbox.pack_start(self.ohdusBox, False, False, 1 )
        for ohdu in self.ohdus:
            toggleButton = Gtk.ToggleButton()
            toggleButton.set_label( '%2s'%(ohdu) )
            toggleButton.set_active(True)
            if ohdu in [11,12,15]: toggleButton.set_active(False)
            #toggleButton.connect('clicked', self.on_ohduButton_clicked )
            self.ohdusBox.pack_start( toggleButton, False, False, 0 )

        refreshButton = Gtk.Button(label='refresh')
        refreshButton.connect('clicked', self.on_ohduButton_clicked )
        subbox.pack_start(refreshButton, False, False, 1)
        return subbox
    
    def build_body(self):
        notebook = Gtk.Notebook()
        notebook.connect( 'switch-page', self.on_switch_page )
        self.currentPageLabel = None
        self.labels = {}
        self.fig = {}
        self.subHeader = {}
        self.imageCanvases = {}
        self.toolbar = {}
        for quantity in self.quantities:
            box = Gtk.VBox()
            #vbox.pack_start( box, True, True, 0)

            self.labels[quantity] = Gtk.Label()
            self.labels[quantity].set_justify(Gtk.Justification.LEFT)
            self.subHeader[quantity] = Gtk.HBox()
            self.resetLabel(quantity)
            self.subHeader[quantity].pack_start(self.labels[quantity], True, True, 0)
            self.fig[quantity] = Figure(dpi=100)
            self.imageCanvases[quantity] = FigureCanvas(self.fig[quantity])
            self.imageCanvases[quantity].set_property('height-request', 500)
            self.fig[quantity].canvas.mpl_connect('button_press_event', self.onclick )
            self.toolbar[quantity] = NavigationToolbar(self.imageCanvases[quantity], self)
            children = self.toolbar[quantity].get_children()
            for i in range(len(children)-3):
                children[i].destroy()
            
            box.pack_start(self.subHeader[quantity], False, False, 0)
            box.pack_start(self.imageCanvases[quantity], True, True, 0)
            box.pack_start(self.toolbar[quantity], False, False, 0)
            notebook.append_page(box,Gtk.Label(label=quantity))
        return notebook
    
    def updateLabel(self, quantity, text ):
        def callback(): 
            self.labels[quantity].set_markup( self.labels[quantity].get_label() + ' ' + text )
            return False
        GLib.idle_add(callback)

    def should_quit(self):
        return self.quit
    
    def resetLabel(self, quantity ):
        def callback():
            percentage = self.percentage_done[quantity]
            if percentage is None:
                self.labels[quantity].set_markup( '<span color="blue">%s</span>'%(quantity) )
                return
            self.labels[quantity].set_markup( '<span color="blue">%s(%d%%)</span>'%(quantity, int(self.percentage_done[quantity]) ) )
            return False
        GLib.idle_add( callback )

    def remove_lock(self, quantity ):
        os.remove( '%s%s.lock'%(self.tablePaths[quantity],self.id) )
        #print 'removing lock', quantity

    def remove_all_locks(self, quantity ):
        locks =  glob.glob('%s*.lock'%(self.tablePaths[quantity]) )
        for lock in locks:
            os.remove( lock )
        print 'all locks removed in thread '+ threading.currentThread().getName() + ' ' + quantity + ' ' + ', '.join(locks)
        self.firstLockMsg[quantity] = True
        
    def create_lock(self, quantity):
        lock_file = '%s%s.lock'%(self.tablePaths[quantity],self.id)
        open( lock_file, 'w' )
        #os.chmod( lock_file, 0o666 )
        return

    def is_locked(self, quantity ):
        lock_files = glob.glob( '%s*.lock'%self.tablePaths[quantity] )
        if len(lock_files) == 0: 
            self.lockTimers[quantity] = None
            self.firstLockMsg[quantity] = True
            return False
        if self.lockTimers[quantity] is None:
            self.lockTimers[quantity] = time.time()
        
        if self.firstLockMsg[quantity]:
            print 'found locks in thread %s %s %s'%(threading.currentThread().getName(), quantity, ', '.join(lock_files))
            self.firstLockMsg[quantity] = False
        pattern = '%s.(.*?)[0-9]*.lock'%self.tablePaths[quantity]
        lock_users = map( lambda lock_file: re.search( pattern, lock_file ).groups()[0], lock_files )
        elapsed_time = min(map( file_age, lock_files ))
        minutes = int(elapsed_time)/60
        hours = minutes/60
        string = '%sm'%minutes if hours == 0 else '%sh%sm'%(hours, minutes)
        self.resetLabel( quantity )
        self.updateLabel(quantity, '<span color="red"><b>is locked by %s for %s</b> will check again soon</span>'%( ','.join(lock_users), string ) )

        if elapsed_time > 5*60 and not self.has_remove_lock_button[quantity]:
            self.has_remove_lock_button[quantity] = True
            print 'should create remove lock button in thread %s %s'%( threading.currentThread().getName(), quantity )
            self.create_remove_lock_button( quantity )
        return True
    
    def start_thread(self, callback, quantity):
        self.thread[quantity] = threading.Thread(target=callback, name=quantity)
        self.thread[quantity].start()
        
    def on_runIDRangeEntry_activate( self, entry ):
        print 'runIDRange activated %s'%entry.get_text()
        #for quantity in self.quantities:
            #self.plot_table( quantity )
        self.plot_table( self.currentPageLabel )
        return
    
    def on_ohduButton_clicked( self, button ):
        #print 'ohdu toggled %s %s'%(button.get_label(), button.get_active())
        button.set_sensitive(False)
        while Gtk.events_pending():
            Gtk.main_iteration()
        #for quantity in self.quantities:
            #self.plot_table( quantity )
        self.plot_table( self.currentPageLabel )
        button.set_sensitive(True)
        #while Gtk.events_pending():
            #Gtk.main_iteration()
        return
        
    def on_omitOutliersButton_clicked( self, button ):
        print 'omitOutliers toggled %s'%button.get_active()
        button.set_sensitive(False)
        while Gtk.events_pending():
            Gtk.main_iteration()
        #for quantity in self.quantities:
            #self.plot_table( quantity )
        self.plot_table( self.currentPageLabel )
        button.set_sensitive(True)
        return
    
    def create_remove_lock_button( self, quantity ):
        def callback():
            print 'removeLock created in thread %s %s'%( threading.currentThread().getName(), quantity )
            removeLockButton = Gtk.Button()
            removeLockButton.set_label('manually remove all %s locks (use with caution)'%(quantity))
            self.subHeader[quantity].pack_start( removeLockButton, False, False, 0 )
            removeLockButton.connect('clicked', self.manually_remove_lock, quantity )
            removeLockButton.show()
            return False
        GLib.idle_add( callback )

    def manually_remove_lock(self, button, quantity ):
        print 'removeLock clicked in thread %s %s'%( threading.currentThread().getName(), quantity )
        shutil.copy2( self.tablePaths[quantity], '%s_%s.backup'%( self.tablePaths[quantity], time.asctime( time.localtime() )) )
        self.remove_all_locks( quantity )
        self.has_remove_lock_button[quantity] = False
        self.resetLabel(quantity)
        GLib.idle_add( button.destroy )
        
    def update_table(self, quantity ):
        '''
        computes an entry for the type specified, reads the previous table and appends to it. During its execution the table file is locked to avoid other instances of the code to modify the file
        '''
        
        if self.is_locked( quantity ): 
            self.wait[quantity] = 5
            return False
        startTime = time.time()
        self.create_lock( quantity )
        
        all_runIDs_reversed = ConniePaths.runID()[::-1]
        if os.path.exists( self.tablePaths[quantity] ):
            data = np.genfromtxt( self.tablePaths[quantity], names=True )
            data = data.tolist()
        else:
            runID = int(getrunIDFromPath(all_runIDs_reversed[0]))
            data = self.make_entry(runID,quantity)
        
        runIDs_done = zip(*data)[0]
        runIDs_done = list(set(runIDs_done))
        self.resetLabel(quantity)
        if len(runIDs_done) == len(all_runIDs_reversed):
            if self.firstAllDoneMsg[quantity]:
                print '%s all done in thread %s'%( threading.currentThread().getName(), quantity )
                self.firstAllDoneMsg[quantity] = False
            self.resetLabel(quantity)
            self.remove_lock(quantity)
            self.updateLabel(quantity, 'all done; sleeping for a minute')
            self.wait[quantity] = 60
            return False
        
        self.firstAllDoneMsg[quantity] = False
        #print len(runIDs_done), len(all_runIDs_reversed)
        self.percentage_done[quantity] = (100*len(runIDs_done))/len(all_runIDs_reversed)
        for runID in all_runIDs_reversed:
            if runID not in runIDs_done:
                print '%s runID %s in thread %s'%(quantity, runID, threading.currentThread().getName() )
                self.resetLabel( quantity )
                self.updateLabel(quantity, 'calculating runID <b>%s</b>'%runID)
                data.extend( self.make_entry(runID, quantity ) )
                break
        
        with open( self.tablePaths[quantity], 'wb' ) as f:
            np.savetxt( f, data, header='#runID ohdu %s'%quantity, fmt='%d %d %.6f' )
        self.remove_lock( quantity )
        self.resetLabel(quantity)
        self.updateLabel(quantity,'concluded <b>%s</b> (%ds)'%(runID, time.time()-startTime) )
        if self.currentPageLabel == quantity:
            self.plot_table( quantity )
        print 'finished %s runID %s in thread %s'%( quantity, runID, threading.currentThread().getName() )
        return False
    
    def make_entry( self, runID, quantity ):
        entry = []
        for ohdu in self.ohdus:
            if self.should_quit():
                self.remove_lock(quantity)
                print 'should quit now from thread %s %s'%( threading.currentThread().getName(), quantity )
                exit(0)
            value = self.compute(runID, ohdu, quantity)
            entry.append( [runID, ohdu, value] )
            self.updateLabel(quantity, str(ohdu))
        return entry
            
    def compute( self, runID, ohdu, quantity ):
        '''
        opens a raw image of the associated runID and ohdu and returns the dark current based on the mean of the diffrence between the medians of the active and overscan area for each line
        '''
        if quantity == 'darkCurrent':
            return ConnieImage.FullImage( runID=runID, ohdu=ohdu, imageType='raw' ).darkCurrentEstimate()
        elif quantity == 'readoutNoise':
            return ConnieImage.FullImage( runID=runID, ohdu=ohdu, imageType='raw' ).readoutNoiseEstimate()
        elif quantity == 'diffMADSqr':
            return ConnieImage.FullImage( runID=runID, ohdu=ohdu, imageType='raw' ).diffMAD()
        else:
            print 'unpredicted quantity', quantity
            exit(1)
    
    def onclick(self, event):
        #print 'onclik signal'
        print 'onclick', (event.x, event.y)
        try:
            ax = event.inaxes
        except:
            return
        lines_points = [ (int(line.get_label()), ax.transData.transform((x,y))) for line in ax.get_lines() for x, y in zip(line.get_xdata(),line.get_ydata()) ]
        #print 'points', points
        picked = [ ( l, ax.transData.inverted().transform((x,y)) ) for l,(x,y) in lines_points if (x-event.x)**2+(y-event.y)**2 < 5**2 ]
        print 'picked', picked
        if len(picked) > 1:
            ohdus = np.unique([ l for l,(x,y) in picked ])
            print 'ohdus', ohdus
            picked_runIDs = [ x for l,(x,y) in picked ]
            runIDmin = min( picked_runIDs )
            runIDmax = max( picked_runIDs )
            self.runIDRangeEntry.set_text( '%d:%d'%( runIDmin, runIDmax ) )
            for button in self.ohdusBox.get_children():
                if not int(button.get_label()) in ohdus: 
                    button.set_active(False)
                else: 
                    button.set_active(True)
            #for quantity in self.quantities:
                #self.plot_table(quantity, max_per_subplot=len(ohdus))
            self.plot_table( self.currentPageLabel, max_per_subplot=len(ohdus))
        elif len(picked) == 1:
            runID = int(picked[0][1][0])
            ohdu = int(picked[0][0])
            print 'should get runID', runID, 'ohdu', ohdu
            open_imageViewer( runID, ohdu )
    
    def on_switch_page( self, notebook, page, page_num ):
        self.currentPageLabel = notebook.get_tab_label_text(page)
        print 'active notebook', self.currentPageLabel
        self.plot_table( self.currentPageLabel )
        #print 'page name', 
        
    def on_mouse_move( self, window, event ):
        #buttonSize = self.omitOutliersButton.get_allocation()
        #print 'button size', buttonSize
        #buttonCoord = window.translate_coordinates( self.omitOutliersButton, event.x, event.y )
        #print 'mouse event omitButton', ( float(buttonCoord[0])/buttonSize[0], float(buttonCoord[1])/buttonSize[1] )
        #print 'mouse event omitButton', ( (event.x-buttonSize.x)/buttonSize.width, (event.y-buttonSize.y)/buttonSize.height )
        canvasSize = self.imageCanvases[self.currentPageLabel].get_allocation()
        print 'canvas size', canvasSize.x, canvasSize.y, canvasSize.width, canvasSize.height
        print 'mouse event fig.canvas', (event.x, event.y)
        print 'axes position', self.grid[0].get_window_extent()
    
    def plot_table( self, quantity, max_per_subplot=4 ):
        '''
        generate image and associate it with the proper imageCanvas object. The range is set by to rangeEntry
        '''
        #if self.plotting[quantity]:
            #self.updateLabel(quantity, '!<span color="yellow">wait for plot</span>' )
            #return
        print 'plot table', quantity
        self.plotting[quantity] = True
        startTime = time.time()
        
        if not os.path.exists( self.tablePaths[quantity] ): return
        data = np.genfromtxt(  self.tablePaths[quantity], names=True )
        
        ohdus = [ int(button.get_label()) for button in self.ohdusBox.get_children() if button.get_active() ]
        runIDRange_str = self.runIDRangeEntry.get_text()
        try:
            runIDMin_str, runIDMax_str = runIDRange_str.split(':')
        except:
            runIDMin_str, runIDMax_str = None, None
        
        runIDMin = utils.try_convertion(runIDMin_str, int)
        runIDMax = utils.try_convertion(runIDMax_str, int)
        
        if runIDMin is None: 
            if self.global_runIDMin is not None:
                runIDMin = min( self.global_runIDMin, data['runID'].min() )
            else:
                runIDMin = data['runID'].min()
            self.global_runIDMin = runIDMin
            self.update = True
        if runIDMax is None: 
            runIDMax = max( self.global_runIDMax, data['runID'].max() )
            self.global_runIDMax = runIDMax
            self.update = True
        

        if runIDMax < 0: runIDMax = data['runID'].max() + runIDMax
        if runIDMin < 0: runIDMin = runIDMax + runIDMin

        if runIDMin > runIDMax:
            print 'wrong zoom %s %s %s'%( quantity, runIDMin, runIDMax )
            return

        runIDMask = np.all( [ data['runID'] >= runIDMin, data['runID'] <= runIDMax], axis=0 )
        
        all_runs = ConniePaths.run()
        all_runBoxes = [ { 'run':run, 'range': (lambda x: [x[0] if x[0] > 1 else x[1], x[-1]])(ConniePaths.runID(run=run)) } for run in all_runs ]
        runBoxes = [ runBox for runBox in all_runBoxes if runBox['range'][1] > runIDMin and runBox['range'][0] < runIDMax ]
        runBoxes[0]['range'][0] = runIDMin
        runBoxes[-1]['range'][1] = runIDMax

        #w, h = self.fig[quantity].get_size_inches()
        #self.fig[quantity].set_size_inches((2*w,h))
        m = int(np.ceil(np.float(len(ohdus))/max_per_subplot))
        #m = int(np.floor( np.sqrt(float(len(ohdus))) ))
        n = float(len(ohdus))/m
        if m==0: return
        for fig in self.fig.values():
            fig.clf()
            fig.set_frameon(True)
        self.grid = []
        ohdus_ranges = []
        for i in range(m):
            share = None if self.grid == [] else self.grid[0]
            #grid.append( self.fig[quantity].add_axes([ .1, .1+(m-i-1)*.8/m, .8, .8/m ], sharex=share, sharey=share ))
            self.grid.append( self.fig[quantity].add_subplot(m,1,i+1, sharex=share, sharey=share ) )
            plt.setp(self.grid[i].get_xticklabels(), visible=False)
            ohdus_ranges.append( ohdus[ int(i*n) : int( (i+1)*n ) ] )

        self.fig[quantity].subplots_adjust(hspace=0)
        ycum = []
        self.lines = {}
        for i, ohdus_range in enumerate(ohdus_ranges):
            for ohdu in ohdus_range:
                for runBox in runBoxes[::-2]:
                    self.grid[i].axvspan( runBox['range'][0], runBox['range'][1], alpha=.02, color='blue' )
                ohduMask = np.all( [runIDMask, data['ohdu']==ohdu], axis=0 )
                y = data[quantity][ohduMask]
                ycum.append(y)
                self.lines[ohdu] = self.grid[i].plot( data['runID'][ohduMask], y, '.', ms=3., label = '%d'%ohdu, rasterized=True )
        
        #print 'lines data', lines[2][0].get_xdata()
        ax = self.grid[0].twiny()
        ax.set_xlim(self.grid[0].get_xlim())
        ax.set_xticks( map( lambda runBox: .5*(runBox['range'][0] + runBox['range'][1]), runBoxes ) )
        ax.set_xticklabels( map( lambda runBox: runBox['run'], runBoxes ) )
        ax.set_zorder(-100)

        medians = np.median(ycum, axis=1)
        mads = stats.MAD(ycum, axis=1)
        if self.omitOutliersButton.get_active():
            val_max = 1.2*medians.max()
        else:
            val_max = 1.05*np.max(ycum)
        
        self.grid[-1].set_xlabel('runID')
        plt.setp(self.grid[-1].get_xticklabels(), visible=True)
        for i in range(m):
            #box = grid[i].get_position()
            self.grid[i].grid(True)
            if quantity == 'darkCurrent':
                self.grid[i].set_ylabel(r'm-m$_{os}$')
            elif quantity == 'readoutNoise':
                self.grid[i].set_ylabel(r'MAD$_{os}$')
            elif quantity == 'diffMADSqr':
                self.grid[i].set_ylabel(r'MAD$^2$-MAD$_{os}^2$')
            self.grid[i].legend( fancybox=True, framealpha=0, bbox_to_anchor=( 1., 1. ), loc='upper left' )
            if m>1:
                self.grid[i].set_ylim((0, val_max))

        #fig.savefig( self.imagePaths[quantity] )
        #plt.close()
        self.updateLabel(quantity, ' <span color="green">plot done (%ds)</span>'%(time.time()-startTime) )

        def callback():
            self.fig[quantity].canvas.draw()
            #self.imageCanvases[quantity].set_from_pixbuf( GdkPixbuf.Pixbuf.new_from_file( self.imagePaths[quantity] ) )
            return False
        
        
        GLib.idle_add(callback)
        self.plotting[quantity] = False
 
