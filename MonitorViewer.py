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

class Logger(object):
    def __init__( self, session, out, key ):
        self.terminal = out
        self.log_file = '/home/mota/public/gui/%s.%s'%(session,key)
        self.key = key
        self.session = session
        
    def write( self, message ):
        self.terminal.write(message)
        log_message = '%s/%s[%s]: %s'%(self.session, self.key, time.strftime("%Y%b%d,%H:%M:%S"), message)
        f = open( self.log_file, 'a+' )
        f.write(message)
        f.close()
        os.chmod( self.log_file, 0o666 )
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
        self.id = '.' + os.environ['HOME'].split('/')[-1]
        self.id += str(int(time.time())%1000000)
        
        print 'id', self.id
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
        for quantity in self.quantities:
            self.plot_table( quantity )

        self.show_all()
        self.quit = False

        for quantity in self.quantities:
            time.sleep(.1)
            self.start_thread( lambda quantity=quantity: self.computation_loop(quantity), quantity )

    def computation_loop( self, quantity ):
        print 'start', quantity, 'in thread loop'
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
            toggleButton.connect('clicked', self.on_ohduButton_clicked )
            self.ohdusBox.pack_start( toggleButton, False, False, 0 )

        return subbox
    
    def build_body(self):
        vbox = Gtk.VBox()
        scroll = Gtk.ScrolledWindow()
        scroll.add_with_viewport(vbox)
        scroll.set_policy(Gtk.PolicyType.NEVER, Gtk.PolicyType.AUTOMATIC)

        self.labels = {}
        self.imageCanvases = {}
        self.subHeader = {}
        for quantity in self.quantities:
            box = Gtk.VBox()
            vbox.pack_start( box, True, True, 0)

            self.labels[quantity] = Gtk.Label()
            self.labels[quantity].set_justify(Gtk.Justification.LEFT)
            self.subHeader[quantity] = Gtk.HBox()
            self.resetLabel(quantity)
            self.subHeader[quantity].pack_start(self.labels[quantity], True, True, 0)
            self.imageCanvases[quantity] = Gtk.Image()
            box.pack_start(self.subHeader[quantity], False, False, 0)
            box.pack_start(self.imageCanvases[quantity], False, False, 0)

        return scroll
        
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
        print 'all locks removed', quantity, ', '.join(locks)
        self.firstLockMsg[quantity] = True
        
    def create_lock(self, quantity):
        lock_file = '%s%s.lock'%(self.tablePaths[quantity],self.id)
        open( lock_file, 'w' )
        os.chmod( lock_file, 0o666 )
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
            print 'found locks', quantity, ', '.join(lock_files)
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
            print 'should create remove lock button', quantity
            self.create_remove_lock_button( quantity )
        return True
    
    def start_thread(self, callback, quantity):
        self.thread[quantity] = threading.Thread(target=callback)
        self.thread[quantity].start()
        
    def on_runIDRangeEntry_activate( self, entry ):
        print 'runIDRange activated', entry.get_text()
        for quantity in self.quantities:
            self.plot_table( quantity )
        return
    
    def on_ohduButton_clicked( self, button ):
        print 'ohdu toggled', button.get_label(), button.get_active()
        button.set_sensitive(False)
        while Gtk.events_pending():
            Gtk.main_iteration()
        for quantity in self.quantities:
            self.plot_table( quantity )
        button.set_sensitive(True)
        #while Gtk.events_pending():
            #Gtk.main_iteration()
        return
        
    def on_omitOutliersButton_clicked( self, button ):
        print 'omitOutliers toggled', button.get_active()
        button.set_sensitive(False)
        while Gtk.events_pending():
            Gtk.main_iteration()
        for quantity in self.quantities:
            self.plot_table( quantity )
        button.set_sensitive(True)
        return
    
    def create_remove_lock_button( self, quantity ):
        def callback():
            print 'removeLock created', quantity
            removeLockButton = Gtk.Button()
            removeLockButton.set_label('manually remove all %s locks (use with caution)'%(quantity))
            self.subHeader[quantity].pack_start( removeLockButton, False, False, 0 )
            removeLockButton.connect('clicked', self.manually_remove_lock, quantity )
            removeLockButton.show()
            return False
        GLib.idle_add( callback )

    def manually_remove_lock(self, button, quantity ):
        print 'removeLock clicked', quantity
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
                print quantity, 'all done'
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
                print quantity, 'runID', runID
                self.resetLabel( quantity )
                self.updateLabel(quantity, 'calculating runID <b>%s</b>'%runID)
                data.extend( self.make_entry(runID, quantity ) )
                break
        
        np.savetxt( self.tablePaths[quantity], data, header='#runID ohdu %s'%quantity, fmt='%d %d %.6f' )
        os.chmod( self.tablePaths[quantity], 0o666 )
        self.remove_lock( quantity )
        self.resetLabel(quantity)
        self.updateLabel(quantity,'concluded <b>%s</b> (%ds)'%(runID, time.time()-startTime) )
        self.plot_table( quantity )
        print 'finished', quantity, 'runID', runID
        return False
    
    def make_entry( self, runID, quantity ):
        entry = []
        for ohdu in self.ohdus:
            if self.should_quit():
                self.remove_lock(quantity)
                print 'should quit now from', quantity
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
    
    def plot_table( self, quantity ):
        '''
        generate image and associate it with the proper imageCanvas object. The range is set by to rangeEntry
        '''
        #if self.plotting[quantity]:
            #self.updateLabel(quantity, '!<span color="yellow">wait for plot</span>' )
            #return
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
        
        if runIDMin < 0: runIDMin = data['runID'].max() + runIDMin
        if runIDMax < 0: runIDMax = data['runID'].max() + runIDMax

        if runIDMin > runIDMax:
            print 'wrong zoom', quantity, runIDMin, runIDMax
            return

        runIDMask = np.all( [ data['runID'] >= runIDMin, data['runID'] <= runIDMax], axis=0 )
        
        all_runs = ConniePaths.run()
        all_runBoxes = [ { 'run':run, 'range': (lambda x: [x[0] if x[0] > 1 else x[1], x[-1]])(ConniePaths.runID(run=run)) } for run in all_runs ]
        runBoxes = [ runBox for runBox in all_runBoxes if runBox['range'][1] > runIDMin and runBox['range'][0] < runIDMax ]
        runBoxes[0]['range'][0] = runIDMin
        runBoxes[-1]['range'][1] = runIDMax

        fig = plt.figure()
        w, h = fig.get_size_inches()
        fig.set_size_inches((2*w,h))
        m = int(np.floor( np.sqrt(float(len(ohdus))) ))
        n = float(len(ohdus))/m
        if m==0: return
        
        grid = []
        ohdus_ranges = []
        for i in range(m):
            share = None if grid == [] else grid[0]
            grid.append( fig.add_axes([ .1, .1+(m-i-1)*.8/m, .8, .8/m ], sharex=share, sharey=share ))
            plt.setp(grid[i].get_xticklabels(), visible=False)
            ohdus_ranges.append( ohdus[ int(i*n) : int( (i+1)*n ) ] )

        ycum = []
        for i, ohdus_range in enumerate(ohdus_ranges):
            for ohdu in ohdus_range:
                for runBox in runBoxes[::-2]:
                    grid[i].axvspan( runBox['range'][0], runBox['range'][1], alpha=.05, color='blue' )
                ohduMask = np.all( [runIDMask, data['ohdu']==ohdu], axis=0 )
                y = data[quantity][ohduMask]
                ycum.append(y)
                grid[i].plot( data['runID'][ohduMask], y, '.', ms=3., label = '%d'%ohdu )
                
        ax = grid[0].twiny()
        ax.set_xlim(grid[0].get_xlim())
        ax.set_xticks( map( lambda runBox: .5*(runBox['range'][0] + runBox['range'][1]), runBoxes ) )
        ax.set_xticklabels( map( lambda runBox: runBox['run'], runBoxes ) )

        medians = np.median(ycum, axis=1)
        mads = stats.MAD(ycum, axis=1)
        if self.omitOutliersButton.get_active():
            val_max = 1.2*medians.max()
        else:
            val_max = 1.05*np.max(ycum)
        
        grid[-1].set_xlabel('runID')
        plt.setp(grid[-1].get_xticklabels(), visible=True)
        for i in range(m):
            #box = grid[i].get_position()
            grid[i].grid(True)
            if quantity == 'darkCurrent':
                grid[i].set_ylabel(r'm-m$_{os}$')
            elif quantity == 'readoutNoise':
                grid[i].set_ylabel(r'MAD$_{os}$')
            elif quantity == 'diffMADSqr':
                grid[i].set_ylabel(r'MAD$^2$-MAD$_{os}^2$')
            grid[i].legend( fancybox=True, framealpha=0, bbox_to_anchor=( 1., 1. ), loc='upper left' )
            grid[i].set_ylim((0, val_max))

        fig.savefig( self.imagePaths[quantity] )
        plt.close()
        self.updateLabel(quantity, ' <span color="green">plot done (%ds)</span>'%(time.time()-startTime) )

        def callback():
            self.imageCanvases[quantity].set_from_pixbuf( GdkPixbuf.Pixbuf.new_from_file( self.imagePaths[quantity] ) )
            return False
        GLib.idle_add(callback)
        self.plotting[quantity] = False
 
