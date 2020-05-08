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

gainAverages = {
    2: 1526,
    3: 1915,
    4: 1637,
    5: 1710,
    6: 1648,
    7: 1498,
    8: 1878,
    9: 1957,
    10: 1810,
    13: 2179,
    14: 1923
    }

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
            #'darkCurrent', 
            #'readoutNoise', 
            #'lambdaBGbin',
            #'lambdaBGbinRAW',
            #'lambdaBGbinBorderRAW',
            'gainElectronbin',
            'gainElectronMAD',
            'gainElectron',
            'gainElectronMAD_line',
            #'diffMedianBG',
            #'diffMADSqr', 
            #'sigmaOS', 
            #'sigmaOSMAD', 
            #'sigmaOSMAD2', 
            #'sigmaOSbin',
            #'gainCu',
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

        #self.plot_table( self.get_active_quantity() )
        self.on_refreshButton_clicked(None)
        
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
        vbox.pack_start(header, expand=False, fill=False, padding=1)
        vbox.pack_start(body, expand=True, fill=True, padding=1)
        return vbox
    
    def build_header(self):
        subbox = Gtk.HBox()

        subbox.pack_start( Gtk.Label( label = 'Dark Current and Readout Noise'), expand=True, fill=True, padding=5 )
        subbox.pack_start( Gtk.Label( label = 'runIDs'), expand=False, fill=False, padding=1 )

        self.runIDRangeEntry = Gtk.Entry()
        if self.range_ is None:
            self.runIDRangeEntry.set_text('-200:auto')
        else:
            self.runIDRangeEntry.set_text(self.range_)
        self.runIDRangeEntry.set_width_chars(len(' auto: auto'))
        subbox.pack_start( self.runIDRangeEntry, expand=False, fill=False, padding=1 )

        ypopover = Gtk.Popover( modal = False, position = Gtk.PositionType.BOTTOM )
        ypopover.add( self.build_yrange_box() )
        self.yrangeButton = Gtk.ToggleButton( label='yrange', relief = Gtk.ReliefStyle.NONE )
        self.yrangeButton.connect('toggled', lambda _: self.show_popover(_, ypopover) )
        
        subbox.pack_start( self.yrangeButton, expand=False, fill=False, padding=1 )
        
        popover = Gtk.Popover( modal = False, position = Gtk.PositionType.BOTTOM )
        popover.add( self.build_ohdu_box() )
        self.ohduButton = Gtk.ToggleButton( label='ohdus', relief = Gtk.ReliefStyle.NONE )
        self.ohduButton.connect('toggled', lambda _: self.show_popover(_, popover) )

        subbox.pack_start(self.ohduButton, expand=False, fill=False, padding=1 )

        qpopover = Gtk.Popover( modal = False, position = Gtk.PositionType.BOTTOM )
        qpopover.add( self.build_quantity_box() )
        self.quantityButton = Gtk.ToggleButton( label='quantity', relief = Gtk.ReliefStyle.NONE )
        self.quantityButton.connect('toggled', lambda _: self.show_popover(_, qpopover) )

        subbox.pack_start( self.quantityButton, expand=False, fill=False, padding=1 )
        
        refreshButton = Gtk.Button( label = 'plot' )
        refreshButton.connect('clicked', self.on_refreshButton_clicked )

        subbox.pack_start(refreshButton, expand=False, fill=False, padding=1)
        return subbox
    
    def show_popover( self, widget, popover, ref=None ):
        if ref is None: ref = widget
        if widget.get_active():
            popover.set_relative_to(ref)
            popover.show_all()
            popover.popup()
        else:
            popover.hide()
        return
    
    def build_explanation_popover( self, widget, ref, quantity ):
        popover = Gtk.Popover( modal = False, position = Gtk.PositionType.LEFT )
        markup = '<span color="blue">%s</span>' % quantity
        label = Gtk.Label()
        popover.add( label )
        label.set_markup( markup )
        widget.connect( 'toggled', lambda _: self.show_popover( widget, popover, ref=ref) )
        return popover

    def build_quantity_box(self):
        #self.quantityBox = 
        vbox = Gtk.VBox()
        self.quantitySelector = []
        for quantity in self.quantities:
            if quantity is self.quantities[0]:
                btn = None
            btn = Gtk.RadioButton.new_with_label_from_widget( btn, quantity )
            self.quantitySelector.append( btn )
            line = Gtk.HBox()
            line.pack_start( btn, expand=False, fill=True, padding=0)
            explanationButton = Gtk.ToggleButton( label = '?', relief = Gtk.ReliefStyle.NONE)
            popover = self.build_explanation_popover( explanationButton, btn, quantity )
            line.pack_end( explanationButton, expand=False, fill=True, padding=0)
            vbox.pack_start( line, expand=False, fill=False, padding=0)
        return vbox
            
    def build_ohdu_box(self):
        self.ohdusBox = Gtk.VBox()
        for ohdu in self.ohdus:
            checkButton = Gtk.CheckButton( label='%02d'%ohdu )
            checkButton.set_active(True)
            if ohdu in [11,12,15]: 
                checkButton.set_active(False)
            self.ohdusBox.pack_start( checkButton, expand=False, fill=False, padding=0 )
        return self.ohdusBox
        
    def build_yrange_box(self):
        self.yrangeBox = Gtk.VBox()
        button = Gtk.RadioButton.new_with_label_from_widget( None, " ymax + 5% " )
        self.yrangeBox.pack_start( button, expand=False, fill=False, padding=0)
        self.yrangeBox.pack_start( 
            Gtk.RadioButton.new_with_label_from_widget( button, " median(y) × 2 " ), expand=False, fill=False, padding=0)
        return self.yrangeBox
        
    def build_body(self):
        box = Gtk.VBox()
        self.imageCanvas = Gtk.Image()
        box.pack_start(self.imageCanvas, expand=True, fill=True, padding=0)
        self.label = Gtk.Label()
        self.subHeader = Gtk.HBox()
        self.subHeader.pack_start( self.label, expand=False, fill=True, padding=0 )
        box.pack_start(self.subHeader, expand=False, fill=False, padding=0)
        return box
    
    def updateLabel(self, quantity, text ):
        def callback(): 
            self.label.set_markup( self.label.get_label() + ' ' + text )
            return False
        GLib.idle_add(callback)

    def should_quit(self):
        return self.quit
    
    def resetLabel(self, quantity ):
        def callback():
            percentage = self.percentage_done[quantity]
            if percentage is None:
                self.label.set_markup( '<span color="blue">%s</span>'%(quantity) )
                return
            self.label.set_markup( '<span color="blue">%s(-%d)</span>'%(quantity, int(self.percentage_done[quantity]) ) )
            return False
        GLib.idle_add( callback )

    def remove_lock(self, quantity ):
        os.remove( '%s%s.lock'%(self.tablePaths[quantity],self.id) )

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
        self.updateLabel(quantity, 
                         '<span color="red"><b>is locked by %s for %s</b> will check again soon</span>'%( ','.join(lock_users), string ) )

        if elapsed_time > 5*60 and not self.has_remove_lock_button[quantity]:
            self.has_remove_lock_button[quantity] = True
            print 'should create remove lock button in thread %s %s'%( threading.currentThread().getName(), quantity )
            self.create_remove_lock_button( quantity, ', '.join(lock_users) )
        return True
    
    def start_thread(self, callback, quantity):
        self.thread[quantity] = threading.Thread(target=callback, name=quantity)
        self.thread[quantity].start()
        
    def on_refreshButton_clicked( self, button ):
        #button.set_sensitive(False)
        self.yrangeButton.set_active(False)
        self.ohduButton.set_active(False)
        self.quantityButton.set_active(False)
        while Gtk.events_pending():
            Gtk.main_iteration()
        self.plot_table( self.get_active_quantity() )
        #button.set_sensitive(True)
        #while Gtk.events_pending():
            #Gtk.main_iteration()
        return
        
    def create_remove_lock_button( self, quantity, users ):
        def callback():
            print 'removeLock created in thread %s %s'%( threading.currentThread().getName(), quantity )
            removeLockButton = Gtk.Button( label='remove %s lock @%s' % (quantity,users) )
            #removeLockButton.set_label('manually remove all %s locks (use with caution)'%(quantity))
            self.subHeader.pack_end( removeLockButton, expand=False, fill=False, padding=0 )
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
            all_runIDs_reversed = ConniePaths.runIDProcessed_v3()[::-1]
        elif quantity in [ 'lambdaBGbin' ]:
            all_runIDs_reversed = self.read_from_table_runIDs( quantity='gainCu' )[::-1]
        else:
            all_runIDs_reversed = ConniePaths.runID()[::-1]

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
                break
        
        if len(data) > 0:
            with open( self.tablePaths[quantity], 'wb' ) as f:
                np.savetxt( f, data, header='#runID ohdu %s'%quantity, fmt='%d %d %.6f' )

        self.remove_lock( quantity )
        self.resetLabel(quantity)
        self.updateLabel(quantity,'concluded <b>%s</b> (%ds)'%(runID, time.time()-startTime) )
        if len(self.eta) > 0:
            estimateEnd = np.mean( self.eta ) * len(runIDs_todo)
            self.updateLabel(quantity,'ETA <b>%s</b>'%( time.ctime( estimateEnd + time.time() ) ) )
        
        if self.get_active_quantity() == quantity: self.plot_table( quantity )
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
        elif quantity == 'gainElectronMAD_line':
            return ConnieImage.FullImage( runID=runID, ohdu=ohdu, imageType='raw', debug=debug ).left().gainElectronMAD_line()
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
            if gainCu is None and 'gainCu' in self.quantities: gainCu = self.read_from_table( runID=runID, ohdu=ohdu, quantity='gainCu' )
            if gainCu is None or gainCu == 0:
                if ohdu in gainAverages.keys():
                    gainCu = gainAverages[ohdu]
                else:
                    return None
            
            gain = gainCu * 3.745e-3 #ADU/keV * keV/e-

            lambdaBGbinRAW = ConnieImage.FullImage( runID=runID, ohdu=ohdu, imageType='raw', debug=debug ).lambdaBGbin( sigma=None, gain=gain, osi=True )
            return lambdaBGbinRAW
        elif quantity == 'lambdaBGbinBorderRAW':
            gainCu = None
            if gainCu is None and 'gainCu' in self.quantities: gainCu = self.read_from_table( runID=runID, ohdu=ohdu, quantity='gainCu' )
            if gainCu is None or gainCu == 0: 
                if ohdu in gainAverages.keys():
                    gainCu = gainAverages[ohdu]
                else:
                    return None
            gain = gainCu * 3.745e-3 #ADU/keV * keV/e-
            
            lambdaBGbinBorderRAW = ConnieImage\
                .FullImage( runID=runID, ohdu=ohdu, imageType='raw', debug=debug )\
                    .left()\
                        .lambdaBGbinBorder( sigma=None, gain=gain )
            return lambdaBGbinBorderRAW
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
        runIDs = list(set([ int(i) for i in data['runID'] ]))
        return runIDs
    
    def get_active_ohdus(self):
        return [ int(button.get_label()) for button in self.ohdusBox.get_children() if button.get_active() ]

    def get_active_yrange(self):
        return [ button.get_label() for button in self.yrangeBox.get_children() if button.get_active() ][0]
    
    def get_active_quantity(self):
        return [ button.get_label() for button in self.quantitySelector if button.get_active() ][0]
    
    def plot_table( self, quantity, max_per_subplot=4 ):
        '''
        generate image and associate it with the proper imageCanvas object. The range is set by to rangeEntry
        '''

        print 'plot table', quantity
        self.plotting[quantity] = True
        startTime = time.time()
        
        self.maximize()
        print 'plotting image'
        allocation = self.imageCanvas.get_allocation()
        print 'allocated size', allocation.width, allocation.height
        
        self.fig = plt.figure( figsize=(allocation.width/100., allocation.height/100.), dpi=100)
        self.imageCanvas.set_from_pixbuf( GdkPixbuf.Pixbuf.new_from_file( '/home/mota/public/gui/loading.gif' ) )
        self.fig.suptitle(quantity)
        
        if not os.path.exists( self.tablePaths[quantity] ): return
        data = np.genfromtxt(  self.tablePaths[quantity], names=True )
        
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

        m = int(np.ceil(np.float(len(ohdus))/max_per_subplot))
        n = float(len(ohdus))/m
        if m==0: return
        self.grid = []
        ohdus_ranges = []
        for i in range(m):
            share = None if self.grid == [] else self.grid[0]
            self.grid.append( self.fig.add_subplot(m,1,i+1, sharex=share, sharey=share ) )
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
                        label = r'${\bf %02d}$ $\mu%.2f$ $\sigma%.2f$' % ( ohdu, np.mean(y), np.std(y) ) )#, rasterized=True )
                self.grid[i].hlines( np.mean(y), min(x), max(x), color=self.lines[ohdu].get_color(), linewidth = .5 )
        
        #print 'lines data', lines[2][0].get_xdata()
        ax = self.grid[0].twiny()
        ax.set_xlim(self.grid[0].get_xlim())
        ax.set_xticks( map( lambda runBox: .5*(runBox['range'][0] + runBox['range'][1]), runBoxes ) )
        ax.set_xticklabels( map( lambda runBox: runBox['run'], runBoxes ) )
        ax.set_zorder(-100)

        val_max = 1
        yrange = self.get_active_yrange()
        if yrange == self.yrangeBox.get_children()[0].get_label():
            val_max = 1.05 * np.nanmax( ycum )
        elif yrange == self.yrangeBox.get_children()[1].get_label():
            medians = np.nanmedian(ycum, axis=1)
            val_max = 2 * np.nanmax( medians )
        
        self.grid[-1].set_xlabel('runID')
        plt.setp(self.grid[-1].get_xticklabels(), visible=True)
        for i in range(m):
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

        try:
            self.fig.tight_layout(rect=(0, 0, .875, .95))
        except:
            print 'could not adjust, please refreh'
        self.fig.subplots_adjust( hspace = 0.05 )
        
        print 'plotting path', self.imagePaths[quantity]
        self.fig.savefig( self.imagePaths[quantity] )
        plt.close()
        
        #self.updateLabel(quantity, ' <span color="green">plot done (%ds)</span>'%(time.time()-startTime) )
        
        def callback():
            print 'plotting image callback'
            self.imageCanvas.set_from_pixbuf( GdkPixbuf.Pixbuf.new_from_file( self.imagePaths[quantity] ) )
            self.imageCanvas.show()
            return False
        
        GLib.idle_add(callback)
        self.plotting[quantity] = False
 
