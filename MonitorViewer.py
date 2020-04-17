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

extraGainRanges = [
    range(12005, 12073),
    range(12073, 12111),
    range(12111, 12213),
    range(12750, 12995),
    ]
extraGainValues = [
    [0, 0, 1535.88, 1925.54, 1635.64, 1715, 1655.42, 1503.19, 1886.29, 1974.87, 1807.92, 0, 0, 2191.9, 1929.86, 2182.48],
    [0, 0, 1533.69, 1924.11, 1641.09, 1717.8, 1648.48, 1508.44, 1885.13, 1973, 1814.79, 0, 0, 2186, 1939.81, 2179.73],
    [0, 0, 1531.06, 1919.26, 1640.19, 1713.7, 1646.31, 1495.45, 1895.85, 1963.93, 1821.96, 0, 0, 2183.45, 1939.44, 2178.98],
    [0, 0, 1530.48, 750.294, 1653.49, 1717.29, 1675.27, 1506.53, 1886.35, 1956.81, 1577.76, 0, 0, 2207.73, 1929.47, 1785.21]
    ]


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
        
    def __init__(self, span=None):
        os.umask(0)

        self.id = '.' + os.environ['HOME'].split('/')[-1]
        self.id += str(int(time.time())%1000000)
        
        print 'id %s'%self.id
        sys.stdout = Logger( self.id, sys.stdout, 'out')
        sys.stderr = Logger( self.id, sys.stderr, 'err')
        
        self.log_file = '/home/mota/public/gui/%s.log'%self.id
        self.quantities = [
            #'lambdaBGbin',
            'lambdaBGbinRAW',
            #'gainElectronbin',
            'gainElectronMAD',
            'gainElectron',
            'diffMedianBG',
            'darkCurrent', 
            'readoutNoise', 
            'diffMADSqr', 
            'sigmaOS', 
            'sigmaOSMAD', 
            'sigmaOSMAD2', 
            'sigmaOSbin',
            'gainCu',
            ]
        
        self.range_ = span
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
        self.show_all()
        self.quit = False
        self.parallel = False

        if self.parallel:
            for quantity in self.quantities:
                print 'request thread for %s %s'%(quantity, threading.active_count())
                self.start_thread( lambda quantity=quantity: self.computation_loop(quantity), quantity )
        else:
            print 'request thread for serial computation %s'%(threading.active_count())
            self.start_thread( self.serial_computation_loop, 'serialComputation' )
        
    def serial_computation_loop( self ):
        self.eta = []
        while 1:
            start = time.time()
            for quantity in self.quantities:
                print 
                print '*** start %s loop in thread %s number of threads %s***'%( quantity, threading.currentThread().getName(), threading.active_count() )
                if self.wait[quantity] > 0:
                    self.wait[quantity] = 0
                    continue
                if self.should_quit():
                    exit(0)
                #try:
                self.update_table( quantity )
                #except:
                    #print "*** exception in update table", quantity, "***"
            self.eta.append( time.time() - start )
        return True

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
        subsubbox = Gtk.HBox()
        scroll = Gtk.ScrolledWindow()
        scroll.set_policy(Gtk.PolicyType.NEVER, Gtk.PolicyType.AUTOMATIC)

        label = Gtk.Label()
        label.set_label('Dark Current and Readout Noise')
        runIDLabel = Gtk.Label()
        runIDLabel.set_label('runIDs')
        self.runIDRangeEntry = Gtk.Entry()
        if self.range_ is None:
            self.runIDRangeEntry.set_text('-100:auto')
        else:
            self.runIDRangeEntry.set_text(self.range_)
        self.runIDRangeEntry.set_width_chars(len(' auto: auto'))

        self.omitOutliersButton = Gtk.ToggleButton()
        self.omitOutliersButton.set_label('outliers')
        self.omitOutliersButton.set_active(False)
        #self.omitOutliersButton.connect( 'clicked', self.on_omitOutliersButton_clicked )
        
        self.interactivePlotButton = Gtk.ToggleButton(label='interactive')
        self.interactivePlotButton.set_active(False)
        #self.interactivePlotButton.connect( 'clicked', self.on_interactivePlotButton_clicked )
        
        subsubbox.pack_start(label, expand=True, fill=True, padding=5 )
        subsubbox.pack_start(runIDLabel, expand=False, fill=False, padding=1 )
        subsubbox.pack_start(self.runIDRangeEntry, expand=False, fill=False, padding=1 )

        #subsubbox.pack_start(yLabel, expand=False, fill=False, padding=1 )
        #subsubbox.pack_start(self.yRangeEntry[, expand=False, fill=False, padding=1 )

        subsubbox.pack_start(self.omitOutliersButton, expand=False, fill=False, padding=1 )
        #subsubbox.pack_start(self.interactivePlotButton, expand=False, fill=False, padding=1 )
        
        popover = Gtk.Popover( modal = False )
        popover.set_position( Gtk.PositionType.BOTTOM )
        popover.add( self.build_ohdu_box() )
        def on_ohduButtonToggled( widget ):
            if widget.get_active():
                popover.set_relative_to(widget)
                popover.show_all()
            else:
                popover.hide()

        ohduButton = Gtk.ToggleButton( label='ohdus' )
        ohduButton.connect('toggled', on_ohduButtonToggled )
        subsubbox.pack_start(ohduButton, expand=False, fill=False, padding=1 )
        
        refreshButton = Gtk.Button(label='refresh')
        refreshButton.connect('clicked', self.on_refreshButton_clicked )
        scroll.add_with_viewport( subsubbox )
        subbox.pack_start(scroll, expand=False, fill=True, padding=1)
        subbox.pack_start(refreshButton, expand=False, fill=False, padding=1)
        return subbox

        
    def build_ohdu_box(self):
        self.ohdusBox = Gtk.VBox()
        for ohdu in self.ohdus:
            checkButton = Gtk.CheckButton( label='%02d'%ohdu )
            checkButton.set_active(True)
            if ohdu in [11,12,15]: 
                checkButton.set_active(False)
            self.ohdusBox.pack_start( checkButton, expand=False, fill=False, padding=0 )
        return self.ohdusBox
        
    
    def build_body(self):
        notebook = Gtk.Notebook()
        notebook.connect( 'switch-page', self.on_switch_page )
        
        notebook.set_scrollable(True)
        notebook.popup_enable()
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
            #self.fig[quantity].canvas.mpl_connect('button_press_event', self.onclick )
            #if self.interactivePlotButton.get_active():
                #self.imageCanvases[quantity] = FigureCanvas(self.fig[quantity])
            #else:
                #self.imageCanvases[quantity] = Gtk.Image()
            self.imageCanvases[quantity] = Gtk.Image()
            
            self.imageCanvases[quantity].set_property('height-request', 500)
            #self.toolbar[quantity] = NavigationToolbar( self.imageCanvases[quantity], self )
            #children = self.toolbar[quantity].get_children()
            #for i in range(len(children)-3):
                #children[i].destroy()
            
            box.pack_start(self.subHeader[quantity], False, False, 0)
            box.pack_start(self.imageCanvases[quantity], True, True, 0)
            #box.pack_start(self.toolbar[quantity], False, False, 0)
            
            notebook.append_page( box, Gtk.Label(label=quantity) )
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
            self.labels[quantity].set_markup( '<span color="blue">%s(-%d)</span>'%(quantity, int(self.percentage_done[quantity]) ) )
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
        
    def on_refreshButton_clicked( self, button ):
        button.set_sensitive(False)
        while Gtk.events_pending():
            Gtk.main_iteration()
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

    def on_interactivePlotButton_clicked( self, button ):
        print 'interactive toggled %s'%button.get_active()
        button.set_sensitive(False)
        while Gtk.events_pending():
            Gtk.main_iteration()
        #self.plot_table( self.currentPageLabel )
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
        
        if quantity in ['sigmaOS', 'sigmaOSbin', 'sigmaOSMAD', 'sigmaOSMAD2', 'gainCu']:
            #all_runIDs_reversed = range(6322, 6521+1)[::-1]
            all_runIDs_reversed = ConniePaths.runIDProcessed_v3()[::-1]
            #print '(update_table)', all_runIDs_reversed
        elif quantity in ['lambdaBGbin', 'lambdaBGbinRAW']:
            #all_runIDs_reversed = extraGainRanges[-1]
            #all_runIDs_reversed.extend( self.read_from_table_runIDs( quantity='gainCu' )[::-1] )
            all_runIDs_reversed = self.read_from_table_runIDs( quantity='gainCu' )[::-1]
        else:
            all_runIDs_reversed = ConniePaths.runID()[::-1]
            #all_runIDs_reversed = range(6322, 6522+1)
        if os.path.exists( self.tablePaths[quantity] ):
            data = np.genfromtxt( self.tablePaths[quantity], names=True )
            data = data.tolist()
        else:
            print '*** path does not exist***'
            lastrunID = all_runIDs_reversed[0]
            data = self.make_entry( lastrunID, quantity )
        
        runIDs_done = []
        if len(data) > 0:
            runIDs_done = zip(*data)[0]
            runIDs_done = list(set(runIDs_done))
            #print '(update_table)runIDs_done=', runIDs_done
            self.resetLabel(quantity)
            runIDs_todo = set(all_runIDs_reversed) - set(runIDs_done)
            if len(runIDs_todo) == 0:
                if self.firstAllDoneMsg[quantity]:
                    print '%s all done in thread %s'%( threading.currentThread().getName(), quantity )
                    self.firstAllDoneMsg[quantity] = False
                self.resetLabel(quantity)
                self.remove_lock(quantity)
                self.updateLabel(quantity, 'all done; sleeping for a minute')
                self.wait[quantity] = 60
                return False
        
        self.firstAllDoneMsg[quantity] = False
        runIDs_todo = set(all_runIDs_reversed) - set(runIDs_done)
        print '(update_table)[%s]len_todo='%quantity, len(runIDs_todo)
        self.percentage_done[quantity] = len(runIDs_todo)
        for runID in list(runIDs_todo)[::-1]:
            #if runID not in runIDs_done:
            print '%s runID %s in thread %s'%(quantity, runID, threading.currentThread().getName() )
            self.resetLabel( quantity )
            self.updateLabel(quantity, 'calculating runID <b>%s</b>'%runID)
            entry = self.make_entry(runID, quantity )
            if len(entry) > 0:
                data.extend( entry )
                break
            else:
                print '(update_table)empty entry at quantity', quantity, 'runID', runID
        
        if len(data) > 0:
            with open( self.tablePaths[quantity], 'wb' ) as f:
                np.savetxt( f, data, header='#runID ohdu %s'%quantity, fmt='%d %d %.6f' )

        self.remove_lock( quantity )
        self.resetLabel(quantity)
        self.updateLabel(quantity,'concluded <b>%s</b> (%ds)'%(runID, time.time()-startTime) )
        if len(self.eta) > 0:
            estimateEnd = np.mean( self.eta ) * len(runIDs_todo)
            self.updateLabel(quantity,'ETA <b>%s</b>'%( time.ctime( estimateEnd + time.time() ) ) )
        
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
            if value is None:
                continue
            entry.append( [runID, ohdu, value] )
            self.updateLabel(quantity, str(ohdu))
        return entry
            
    def compute( self, runID, ohdu, quantity, debug=False ):
        '''
        opens a raw image of the associated runID and ohdu and returns the dark current based on the mean of the diffrence between the medians of the active and overscan area for each line
        '''
        if debug:
            print '(MonitorViewer.compute)quantity=%s'%quantity
        if quantity == 'darkCurrent':
            return ConnieImage.FullImage( runID=runID, ohdu=ohdu, imageType='raw' ).darkCurrentEstimate()
        elif quantity == 'readoutNoise':
            return ConnieImage.FullImage( runID=runID, ohdu=ohdu, imageType='raw' ).readoutNoiseEstimate()
        elif quantity == 'diffMADSqr':
            return ConnieImage.FullImage( runID=runID, ohdu=ohdu, imageType='raw' ).diffMAD()
        elif quantity == 'diffMedianBG':
            return ConnieImage.FullImage( runID=runID, ohdu=ohdu, imageType='raw', debug=debug ).diffMedianBG()
        elif quantity == 'gainElectron':
            return ConnieImage.FullImage( runID=runID, ohdu=ohdu, imageType='raw', debug=debug ).left().gainElectron()
        elif quantity == 'gainElectronbin':
            return ConnieImage.FullImage( runID=runID, ohdu=ohdu, imageType='raw', debug=debug ).left().gainElectronbin()
        elif quantity == 'gainElectronMAD':
            return ConnieImage.FullImage( runID=runID, ohdu=ohdu, imageType='raw', debug=debug ).left().gainElectronMAD()
        elif quantity == 'sigmaOS':
            return ConnieImage.FullImage( runID=runID, ohdu=ohdu, imageType='scn', debug=debug ).sigmaOS()
        elif quantity == 'sigmaOSbin':
            return ConnieImage.FullImage( runID=runID, ohdu=ohdu, imageType='scn', debug=debug ).sigmaOSbin()
        elif quantity == 'sigmaOSMAD':
            return ConnieImage.FullImage( runID=runID, ohdu=ohdu, imageType='scn', debug=debug ).sigmaOSMAD()
        elif quantity == 'sigmaOSMAD2':
            return ConnieImage.FullImage( runID=runID, ohdu=ohdu, imageType='scn', debug=debug ).sigmaOSMAD2()
        elif quantity == 'gainCu':
            #print '(compute)runID=', runID, 'ohdu=', ohdu, 'quantity=', quantity
            gainCu = ConnieImage.read_gainCu( runID=runID, ohdu=ohdu )
            #print '(compute)gainCu=', gainCu
            return gainCu
        elif quantity == 'lambdaBGbin':
            print 'lambdaBGbin', 'runID=', runID, 'ohud=', ohdu
            if runID == 9777:
                return None
            gainCu = None
            for index, extraGainRange in enumerate(extraGainRanges):
                if runID in extraGainRange:
                    gainCu = extraGainValues[index][ohdu]
                    if gainCu is None: return None
                    if gainCu == 0: return None
                    break
            if gainCu is None and 'gainCu' in self.quantities: gainCu = self.read_from_table( runID=runID, ohdu=ohdu, quantity='gainCu' )
            if gainCu is None or gainCu == 0: return None
            gain = gainCu * 3.745e-3 #ADU/keV * keV/e-
            
            #try:
            lambdaBGbin = ConnieImage.FullImage( runID=runID, ohdu=ohdu, imageType='scn', debug=debug ).lambdaBGbin( sigma=None, gain=gain )
            #except:
                #print "*** lambdaBGbin error for", runID, ohdu, "***"
                #return None
            return lambdaBGbin
        elif quantity == 'lambdaBGbinRAW':
            print 'lambdaBGbin', 'runID=', runID, 'ohud=', ohdu
            gainCu = None
            #for index, extraGainRange in enumerate(extraGainRanges):
                #if runID in extraGainRange:
                    #gainCu = extraGainValues[index][ohdu]
                    #if gainCu is None: return None
                    #if gainCu == 0: return None
                    #break
            if gainCu is None and 'gainCu' in self.quantities: gainCu = self.read_from_table( runID=runID, ohdu=ohdu, quantity='gainCu' )
            if gainCu is None or gainCu == 0: return None
            gain = gainCu * 3.745e-3 #ADU/keV * keV/e-

            lambdaBGbinRAW = ConnieImage.FullImage( runID=runID, ohdu=ohdu, imageType='raw', debug=debug ).lambdaBGbin( sigma=None, gain=gain, osi=True )
            return lambdaBGbinRAW
        else:
            print 'unpredicted quantity', quantity
            exit(1)
    
    def read_from_table( self, runID, ohdu, quantity ):
        if os.path.exists( self.tablePaths[quantity] ):
            data = np.genfromtxt( self.tablePaths[quantity], names=True )
        else:
            return None
        #print '(MonitorViewer.read_from_table)quantity=', quantity
        #print '(MonitorViewer.read_from_table)runID=', runID, 'ohdu=', ohdu
        #data = data[ data['runID'] == runID ]
        #print '(MonitorViewer.read_from_table)data=', data
        ret = data[ np.all([data['runID']==runID, data['ohdu']==ohdu], axis=0) ]
        if len(ret) == 0: return None
        print '(MonitorViewer.read_from_table)ret=', ret[0][-1]
        return ret[0][-1]

    def read_from_table_runIDs( self, quantity ):
        if os.path.exists( self.tablePaths[quantity] ):
            data = np.genfromtxt( self.tablePaths[quantity], names=True )
        else:
            return None
        #print '(update_table)', data['runID']
        runIDs = list(set([ int(i) for i in data['runID'] ]))
        #print '(update_table)', runIDs
        return runIDs
    
    #def onclick(self, event):
        ##print 'onclik signal'
        #print 'onclick', (event.x, event.y)
        #try:
            #ax = event.inaxes
        #except:
            #return
        #lines_points = [ (int(line.get_label()), ax.transData.transform((x,y))) for line in ax.get_lines() for x, y in zip(line.get_xdata(),line.get_ydata()) ]
        ##print 'points', points
        #picked = [ ( l, ax.transData.inverted().transform((x,y)) ) for l,(x,y) in lines_points if (x-event.x)**2+(y-event.y)**2 < 5**2 ]
        #print 'picked', picked
        #if len(picked) > 1:
            #ohdus = np.unique([ l for l,(x,y) in picked ])
            #print 'ohdus', ohdus
            #picked_runIDs = [ x for l,(x,y) in picked ]
            #runIDmin = min( picked_runIDs )
            #runIDmax = max( picked_runIDs )
            #self.runIDRangeEntry.set_text( '%d:%d'%( runIDmin, runIDmax ) )
            #for button in self.ohdusBox.get_children():
                #if not int(button.get_label()) in ohdus: 
                    #button.set_active(False)
                #else: 
                    #button.set_active(True)
            #self.plot_table( self.currentPageLabel, max_per_subplot=len(ohdus))
        #elif len(picked) == 1:
            #runID = int(picked[0][1][0])
            #ohdu = int(picked[0][0])
            #print 'should get runID', runID, 'ohdu', ohdu
            #open_imageViewer( runID, ohdu )
    
    def on_switch_page( self, notebook, page, page_num ):
        self.currentPageLabel = notebook.get_tab_label_text(page)
        print
        print 'active notebook', self.currentPageLabel
        self.plot_table( self.currentPageLabel )
        
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
    
    def get_active_ohdus(self):
        return [ int(button.get_label()) for button in self.ohdusBox.get_children() if button.get_active() ]
    
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
        
        if self.interactivePlotButton.get_active():
            print 'plotting interactive'
            self.fig[quantity] = Figure(dpi=100)
            self.imageCanvases[quantity] = FigureCanvas(self.fig[quantity])
            
            #self.toolbar[quantity] = NavigationToolbar(self.imageCanvases[quantity], self)
            #children = self.toolbar[quantity].get_children()
            #for i in range(len(children)-3):
                #children[i].destroy()
        else:
            print 'plotting image'
            allocation = self.imageCanvases[quantity].get_allocation()
            print 'allocated size', allocation.width, allocation.height
            
            self.fig[quantity] = plt.figure( figsize=(allocation.width/100., allocation.height/100.), dpi=100)
            self.imageCanvases[quantity].set_from_pixbuf( GdkPixbuf.Pixbuf.new_from_file( '/home/mota/public/gui/loading.gif' ) )
            #self.imageCanvases[quantity] = Gtk.Image()

            #self.toolbar[quantity] = None
        
        if not os.path.exists( self.tablePaths[quantity] ): return
        data = np.genfromtxt(  self.tablePaths[quantity], names=True )
        
        #ohdus = [ int(button.get_label()) for button in self.ohdusBox.get_children() if button.get_active() ]
        ohdus = self.get_active_ohdus()
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
            #runIDMax = max( self.global_runIDMax, data['runID'].max() )
            runIDMax = data['runID'].max()
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
        #for fig in self.fig[quantity].values():
            #fig.clf()
            #fig.set_frameon(True)
        self.grid = []
        ohdus_ranges = []
        for i in range(m):
            share = None if self.grid == [] else self.grid[0]
            #grid.append( self.fig[quantity].add_axes([ .1, .1+(m-i-1)*.8/m, .8, .8/m ], sharex=share, sharey=share ))
            self.grid.append( self.fig[quantity].add_subplot(m,1,i+1, sharex=share, sharey=share ) )
            plt.setp( self.grid[i].get_xticklabels(), visible=False )
            ohdus_ranges.append( ohdus[ int(i*n) : int( (i+1)*n ) ] )

        ycum = []
        self.lines = {}
        for i, ohdus_range in enumerate(ohdus_ranges):
            for ohdu in ohdus_range:
                for runBox in runBoxes[::-2]:
                    self.grid[i].axvspan( runBox['range'][0], runBox['range'][1], alpha=.02, color='blue' )
                ohduMask = np.all( [runIDMask, data['ohdu']==ohdu], axis=0 )
                y = data[quantity][ohduMask]
                if len(y) == 0:
                    continue
                ycum.append(y)
                x = data['runID'][ohduMask]
                self.lines[ohdu], = self.grid[i].plot( x, y, '.', ms=3., 
                        label = '${\bf %02d}$ $\mu%.2f$ $\sigma%.2f$' % ( ohdu, np.mean(y), np.std(y) ) )#, rasterized=True )
                self.grid[i].hlines( np.mean(y), min(x), max(x), color=self.lines[ohdu].get_color() )
        
        #print 'lines data', lines[2][0].get_xdata()
        ax = self.grid[0].twiny()
        ax.set_xlim(self.grid[0].get_xlim())
        ax.set_xticks( map( lambda runBox: .5*(runBox['range'][0] + runBox['range'][1]), runBoxes ) )
        ax.set_xticklabels( map( lambda runBox: runBox['run'], runBoxes ) )
        ax.set_zorder(-100)

        #mads = stats.MAD(ycum, axis=1)

        print 'ycum length', len(ycum)
        val_max = 1
        if len(ycum) > 1:
            if len(ycum[0]) > 0: 
                if self.omitOutliersButton.get_active():
                    medians = np.nanmedian(ycum, axis=1)
                    val_max = 2 * np.nanmax( medians )
                else:
                    val_max = 1.05 * np.nanmax( ycum )
                #val_max = 1.05*np.nanmax(ycum)
        
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
                self.grid[i].set_ylim(( (1-.05)*np.nanmin(ycum), val_max ))

        self.fig[quantity].tight_layout(rect=(0, 0, .875, 1))
        self.fig[quantity].subplots_adjust( hspace = 0.05 )
        
        if not self.interactivePlotButton.get_active():
            print 'plotting path', self.imagePaths[quantity]
            self.fig[quantity].savefig( self.imagePaths[quantity] )
            plt.close()
        
        #self.updateLabel(quantity, ' <span color="green">plot done (%ds)</span>'%(time.time()-startTime) )
        
        def callback():
            if self.interactivePlotButton.get_active():
                print 'plotting interactive callback'
                self.fig[quantity].canvas.draw()
            else:
                print 'plotting image callback'
                self.imageCanvases[quantity].set_from_pixbuf( GdkPixbuf.Pixbuf.new_from_file( self.imagePaths[quantity] ) )
                self.imageCanvases[quantity].show()
            return False
        
        GLib.idle_add(callback)
        self.plotting[quantity] = False
 
