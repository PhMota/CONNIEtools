#!/usr/bin/env python
import ROOT
from ROOT import AddressOf
#import ctypes
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.colors
import matplotlib.patches as patches
import numpy as np
import scipy.weave
#import array
#import copy
from scipy import stats
import time, sys, os, re
import operator
from timer import timer
import root2numpy
import root_numpy
#import rootpy
import locker
import shutil
import json
import io
import codecs
import ast

#plt.rc('text', usetex=True)

#def update_mean_std( mean, sigma, N, entry ):
    #updated_mean = (mean*N + entry)/float(N+1)
    #updated_sigma = np.sqrt( ( N*( sigma**2 + mean**2 ) + entry**2 )/float(N+1) - updated_mean**2 )
    #return updated_mean, updated_sigma

iaddattr = lambda object, name, value: setattr( object, name, getattr( object, name ) + value )

def std( x, weights ):
    weights[weights!=weights] = 0
    std2 = np.average( x**2, weights = weights ) - np.average( x, weights = weights )**2
    return np.sqrt( np.abs(std2) )

def cov2(m, a):
    ddof = 0# len(a)-1
    v1 = np.sum(a) # N
    v2 = np.sum(a**2) # N
    m -= np.sum(m * a, axis=1, keepdims=True) / v1 # m = m - (sum m)/N
    return np.dot(m * a, m.T) * v1 / (v1**2 - ddof * v2) # m*m.T * N/(N**2 - (N-1)*N)

def cov(m, a):
    v1 = np.sum(a) # N
    m -= np.sum(m * a, axis=1, keepdims=True) / v1 # m = m - (sum m)/N
    return np.dot(m * a, m.T) / v1 # m*m.T * N/(N**2 - (N-1)*N)

def gaussian( point, center, sigmas, directions = None, integral = None, amplitude = None ):
    if directions != None:
        #projected_dist = project( point - center, v = directions )/sigmas**2
        projected_dist = np.inner( directions[:,0], np.array(point) - center )**2/sigmas[0]**2 + np.inner( directions[:,1], np.array(point) - center )**2/sigmas[1]**2
    else:
        projected_dist = np.sum( np.array(point) - center )**2/sigmas[0]**2
    if integral != None:
        return integral/(2.*np.pi*sigmas[0]*sigmas[1]) * np.exp( -.5 * projected_dist)
    elif amplitude != None:
        #print "gaussian", amplitude, np.exp( -.5 * projected_dist)
        return amplitude * np.exp( -.5 * projected_dist)
    print "must give either integral or amplitude"
    return 0

def project( x, y = None, v = None ):
    if y != None:
        x = np.array([x,y])
    return np.dot( v.T, x )

def make_gaussian2D( x, y, weights = None ):
    """
    return a gaussian function based on the input 2D distribution, the center, the covariance and the directions
    """
    if weights != None:
        ddof = len(weights)-1
        total_weight = np.sum(weights)
        total_weight_sqr = np.sum(weights**2)
        average = np.average([x,y], weights=weights, axis = 1)
        shifted_xy = [x,y] - average[:,None]
        cov = np.dot( shifted_xy * weights, shifted_xy.T )*total_weight/(total_weight**2 - ddof*total_weight_sqr )
        integral = total_weight
    else:
        average = np.average([x,y], axis = 1)
        cov = np.cov( [x,y], shifted_xy.T )
        integral = len(x)
    
    #average = np.average([x,y], weights=weights, axis = 1)
    #shifted_xy = [x,y] - average[:,None]
    #cov = np.dot( shifted_xy * weights, shifted_xy.T )*total_weight/(total_weight**2 - ddof*total_weight_sqr )
    eigenvalues, eigenvectors = np.linalg.eig( cov )
    #print shifted_xy.shape, eigenvectors.shape
    projected_xy = project( shifted_xy, v = eigenvectors )
    variance = np.average( projected_xy**2, weights=weights, axis=1)
    #print "projected", projected_xy[:,0]
    #print "innerSingle", np.inner( eigenvectors[:,0], np.array([x[0],y[0]]) - average ), np.inner( eigenvectors[:,1], np.array([x[0],y[0]]) - average )
    #print "var", np.sqrt(variance), np.sqrt(np.abs(eigenvalues)), std( projected_xy[0,:], e ), np.std( projected_xy[0,:] )
    return lambda xy: gaussian(xy, average, np.sqrt(variance), eigenvectors, integral = integral ), average, np.sqrt(variance), eigenvectors

def parse_tests(text):
    """
	reads:
	"muon(r): $stdT < .5 and $nSavedPix > 200
	electron(k): $rmseT > 50 and $stdT > .8 ]"
	and parses to:
	[ ["muon", ['r', "$stdT < .5 and $nSavedPix > 200" ]],
	  ["electron", ['k', "$rmseT > 50 and $stdT > .8" ]] ]
    """
    lines = text.split('\n')
    print "lines", lines
    parsed_lines = []
    for line in lines:
	m = re.match( "^(.+?)\((\w+)\):\ (.*)$", line )
	print [ [m.group(1), [m.group(2), m.group(3) ]] ]
	parsed_lines += [ [m.group(1), [m.group(2), m.group(3) ]] ]
    return parsed_lines

def string_into_test( obj, s, list_of_vars ):
    s_new = s
    #for var in sorted(self.vars_scalar + self.vars_extra, key = len, reverse=True ):
    for var in list_of_vars:
	s_new = re.sub("\$%s"%var, "%s.%s"%(obj,var), s_new)
    if s_new == s: return None
    if "$" in s_new: return None
    return s_new

def string_into_test2( obj, s, list_of_vars ):
    s_new = s
    for var in list_of_vars:
	s_new = re.sub("\$%s"%var, "%s.%s[n]"%(obj,var), s_new)
    if s_new == s: return None
    if "$" in s_new: return None
    return s_new

varsCor = ['logEloss', 'logRedRmseT', 'redn0', 'relLen', 'redVarT', 'logRedVarE' ]

def updateTrainFile( username, base_catalog, train_catalog, track_id, category, guess = None ):
    t = timer('updateTrainingFile')
    train_file = base_catalog.split('.root')[0] + '.learn.' + username + '.dict'
    if not os.path.exists( train_file ):
        print 'some error', train_file, 'does not exist'
        return 'doesnotexists'
    print username, track_id, category, guess
    
    track_id = int(track_id)
    if not os.path.exists( train_file ):
        print 'file does not exist', train_file
        exit(0)
    categories, done, data, data2, properties = read_categories( username, base_catalog, train_catalog )
    if not category in data:
        data[category] = []
    # if guess is given, check for the actual guess
    guess = eval(guess)
    if not guess is None:
        # if guess is not the category, update with a miss in the guessed category
        if guess[0] != category:
            entry = { 'missed': track_id, 'train': category }
            data[ guess[0] ] += [ entry ]
    #update the category with the track_id and the updated properties of this category
    #last_entry = data[ category ][-1] if len( data[category] ) > 1 else None
    #properties = None
    this_track = root_numpy.root2array( base_catalog, treename = 'hitSumm', branches = varsCor, start = track_id-1, stop = track_id )
    
    try:
        prop = dict(properties[category])
        N = prop['N']
        prop['N'] += 1
        for i, var in enumerate(varsCor):
            prop['<%s>'%var] = (prop['<%s>'%var]*N + float(this_track[var][0]))/(N+1)
            prop['<%s.2>'%var] = (prop['<%s.2>'%var]*N + float(this_track[var][0]**2))/(N+1)
            for var2 in varsCor[:i]:
                prop['<%s.%s>'%(var,var2)] = (prop['<%s.%s>'%(var,var2)]*N + float(this_track[var][0]*this_track[var2][0]))/(N+1)
    except:
        properties[category] = {}
        prop = properties[category]
        prop = {}
        N = 1
        prop['N'] = 1
        for i, var in enumerate(varsCor):
            prop['<%s>'%var] = float(this_track[var][0])
            prop['<%s.2>'%var] = float(this_track[var][0]**2)
            for var2 in varsCor[:i]:
                prop['<%s.%s>'%(var,var2)] = float(this_track[var][0]*this_track[var2][0])
    
    entry = {'trackID': track_id, 'guess': guess, 'properties': prop }
    if not track_id in data[category] and not entry in data[category]:
        data[category] += [ entry ]

    shutil.copy( train_file, train_file+'~' )
    sdata = json.dumps( data, indent = 2, sort_keys=True, ensure_ascii=False )
    with io.open( train_file, 'w', encoding='utf8') as json_file:
        json_file.write(unicode(sdata))
    return None


def read_categories( username, base_catalog, train_catalog, keys = None ):
    t = timer('read_categories')
    numberOfTrained = 0
    done = []
    #todo = [ True for i in range(len(self.nSavedPix)) ]
    #filename = self.fname.split('.root')[0] + '.learn.%s.dict'%user
    
    categories = {}
    total = {}
    guess = {}
    correct = {}
    efficiency = {}
    contamination = {}
    properties = {}
    falsePositive = {}
    falseNegative = {}
    cumFalsePositive = {}
    cumFalseNegative = {}
    lastFalsePositive = {}
    lastFalseNegative = {}
    matrixFalsePositive = {}
    matrixFalseNegative = {}
    
    lastN = 100
    train_file = base_catalog.split('.root')[0] + '.learn.' + username + '.dict'
    if os.path.exists(train_file):
        with codecs.open(train_file, 'rU', 'utf-8') as data_file:
            try:
                data = json.load( data_file, encoding = 'utf-8' )
            except:
                data = {}
        
        for key, values in data.iteritems():
            if not keys is None: 
                if not key in keys: continue
            if not key in categories: properties[key] = {}
            if not key in categories: categories[key] = []
            if not key in total: total[key] = 0
            if not key in guess: guess[key] = 0
            if not key in correct: correct[key] = 0
            falsePositive[key] = []
            falseNegative[key] = []
            
            
            for value in values:
                try:
                    if 'missed' in value:
                        guess[key] += 1
                        falsePositive[key] += [ value['train'] ]
                        continue
                except: 
                    pass
                falsePositive[key] += [0.]
                total[key] += 1
                try: 
                    trackID = value['trackID']
                    try:
                        if value['guess'][0] == key: 
                            guess[key] += 1
                            correct[key] += 1
                            falseNegative[key] += [0.]
                        else:
                            falseNegative[key] += [ value['guess'][0] ]
                    except:
                        pass
                except: 
                    trackID = value
                
                
                done += [trackID]
                numberOfTrained += 1
                
                categories[key] += [ trackID ]
                try:
                    properties[key] = value['properties']
                except:
                    properties[key] = value['category_properties']
        #print 'efficiency and contamination'
        
        for key, values in data.iteritems():
            matrixFalsePositive[key] = { key2: 
                                    np.mean( map( bool, falsePositive[key][ max(len(falsePositive[key])-lastN, 0): ]) ) if key2 is key else
                                    np.mean( map( lambda s: 1. if s is key2 else 0, falsePositive[key][ max(len(falsePositive[key])-lastN, 0): ]) )
                                    for key2 in data.keys() }
            matrixFalseNegative[key] = { key2: 
                                    np.mean( map( bool, falseNegative[key][ max(len(falseNegative[key])-lastN, 0): ]) ) if key2 is key else
                                    np.mean( map( lambda s: 1. if s is key2 else 0, falseNegative[key][ max(len(falseNegative[key])-lastN, 0): ]) )
                                    for key2 in data.keys() }

            #print properties[key]
            #print key, guess[key], total[key], correct[key]
            cumFalseNegative[key] = np.mean( map( bool, falseNegative[key] ) )
            cumFalsePositive[key] = np.mean( map(bool, falsePositive[key] ) )
            lastFalsePositive[key] = np.mean( map( bool, falsePositive[key][ max(len(falsePositive[key])-lastN, 0): ]) )
            lastFalseNegative[key] = np.mean( map( bool, falseNegative[key][ max(len(falseNegative[key])-lastN, 0): ]) )
            efficiency[key] = 1
            contamination[key] = 0
            if guess[key] > 0 and total[key] > 0:
                #print 'efficiency', 1. - (guess[key] - correct[key])/float(guess[key]), 'contamination', (total[key] - correct[key])/float(total[key])
                #efficiency[key] = 1. - (guess[key] - correct[key])/float(guess[key])
                #contamination[key] = (total[key] - correct[key])/float(total[key])
                efficiency[key] = 1. - (total[key] - correct[key])/float(total[key])
                contamination[key] = (guess[key] - correct[key])/float(guess[key])

    else:
        open(train_file,'a').close()
        data = {}
        numberOfTrained = 0

    #print "number of trained", numberOfTrained

    print matrixFalseNegative
    data2 = {}
    data2['total'] = total
    data2['guess'] = guess
    data2['correct'] = correct
    data2['efficiency'] = efficiency
    data2['contamination'] = contamination
    data2['cumFalseNegative'] = cumFalseNegative
    data2['cumFalsePositive'] = cumFalsePositive
    data2['lastFalseNegative'] = lastFalseNegative
    data2['lastFalsePositive'] = lastFalsePositive
    data2['matrixFalseNegative'] = matrixFalseNegative
    data2['matrixFalsePositive'] = matrixFalsePositive
    #print 'data', data
    return categories, done, data, data2, properties


####################################################################################################################
####################################################################################################################
####################################################################################################################
class Catalog2:
####################################################################################################################
    def listRunID( self, fname, nEntries ):
        print "opening", fname
        tfile = ROOT.TFile.Open( fname )
        config = tfile.Get("config")
        tree = tfile.Get("hitSumm")
        print tree
        ROOT.gROOT.ProcessLine("\
            struct event_t {\
            Int_t runID;\
            }" )
        event = ROOT.event_t()
        tree.SetBranchAddress( "runID", AddressOf( event, "runID" ) )

        runID = {}
        for i in xrange( tree.GetEntries() ):
            tree.GetEntry(i)
            if not event.runID in runID: runID[event.runID] = 0
            else: runID[event.runID] += 1

	tree.ResetBranchAddress(0)
        #print runID.items()
        print "runID", min(runID.keys()), max(runID.keys()), len(runID.keys())

        #nEntries = 50
        step = 0
        minIndex = 0
        while minIndex < max(runID.keys()):
            cutfname = fname.split('.root')[0] + 'cut%03i.root' % step
            cuttfile = ROOT.TFile.Open( cutfname, "RECREATE" )
            config.CopyTree("")
            minIndex = min(runID) + step*nEntries
            maxIndex = min(runID) + (step+1)*nEntries
            #tree.CopyTree("flag == 0 && nSavedPix > 0 && runID >= %i && runID < %i" % ( minIndex, maxIndex ) )
            tree.CopyTree("runID >= %i && runID < %i" % ( minIndex, maxIndex ) )
            cuttfile.Write()
            cuttfile.Close()
            step += 1
            minIndex = min(runID) + step*nEntries

	#tree.Delete()
        tfile.Close()
        print step
        return step
    
####################################################################################################################
    def splitCatalog( self, fname ):
        print "opening", fname
        t = timer('split')
        
        data = root_numpy.root2array(fname, treename = 'hitSumm', branches = ['runID', 'ohdu'] ).view(np.recarray)
        
        runIDs = set(data['runID'])
        
        tfile = ROOT.TFile.Open( fname )
        config = tfile.Get("config")
        tree = tfile.Get("hitSumm")
        print tree
        
        for runID in runIDs:
            t = timer('runID')
            print runID
            cutfname = fname.split('.root')[0] + 'runID%04i.root' % runID
            if os.path.exists(cutfname):
                print cutfname, 'already exists; not overwriting'
                continue
            cuttfile = ROOT.TFile.Open( cutfname, "RECREATE" )
            config.CopyTree("")
            tree.CopyTree("runID == %i" % runID )
            cuttfile.Write( '', tree.kOverwrite )
            cuttfile.Close()
            self.addID( cutfname )
        tfile.Close()
        del t
        return
####################################################################################################################
    def addID( self, fname ):
        print "opening", fname
        if 'trackID' in root_numpy.list_branches( fname, treename = 'hitSumm' ):
            tfile = ROOT.TFile.Open( fname )
            tree = tfile.Get('hitSumm')
            print tree.Print()
            print 'trackID is already at', fname
            data = root_numpy.root2array(fname, treename = 'hitSumm', branches = ['trackID', 'eventID'] ).view(np.recarray)
            print data['trackID']
            print data['eventID']
            return 0
        #t = timer('split')

        data = root_numpy.root2array(fname, treename = 'hitSumm', branches = ['runID', 'ohdu'] ).view(np.recarray)
        runIDs = set(data['runID'])
        
        trackID = np.array( range( 1, len(data['runID']) + 1 ) )
        eventID = np.array( [ '%04d%02d%09d'%(r,o,t) for r, o, t in zip( data['runID'], data['ohdu'], trackID ) ] )
        trackID = trackID.view(np.recarray)
        eventID = eventID.view(np.recarray)
        arr = np.array( zip(trackID, eventID), dtype=[('trackID', np.int32), ('eventID', 'S32')] ).view(np.recarray)
        
        root_numpy.array2root( arr, fname, treename = 'hitSumm', mode='update' )
        
        return
####################################################################################################################
####################################################################################################################
####################################################################################################################
#class Catalog3:
    #def updateCatalog( self, fname ):
        #print "opening", fname
        #t = timer('update')
        
        #data = root_numpy.root2array(fname, treename = 'hitSumm', branches = ['runID', 'ohdu'] ).view(np.recarray)
####################################################################################################################
####################################################################################################################
####################################################################################################################
class Catalog:
####################################################################################################################
    def readTreeFromFile( self, fname ):
        print 'reading catalog', fname
        #return ROOT.TFile.Open('/home/mota/Projects/Connie/damic/damicViewer/catalog_55Fe.root').Get("hitSumm")
        return ROOT.TFile.Open( fname, 'update' ).Get("hitSumm")
####################################################################################################################
    def __init__( self, fname = None ):
        t = timer('Catalog.__init__')
        self.DATATREE = 'hitSumm'
        self.CONFTREE = 'config'

        self.vars_scalar = ['flag', 'runID','nSavedPix','ohdu','E0','n0']
        self.vars_list = ['xPix','yPix','ePix','level']
        self.vars_extra = [
            'L', 'T',
            'stdL', 'stdT',
            'rmseL', 'rmseT',
                ]
        #self.varsCor = ['logEloss', 'logRedRmseT', 'redn0', 'relLen', 'redVarT', 'logRedVarE' ]
        
        self.vars_all = self.vars_scalar + self.vars_list
        self.vars_listextra = self.vars_all + self.vars_extra

        self.fname = fname
        if fname == None:
            return
        print fname
        if not os.path.exists( fname ):
            print "catalog not found", fname
        #f = ROOT.TFile( fname )
        #self.tree = f.Get("hitSumm")
        #if not complete:
            #self.parseTree( self.tree )
        #else:
            #print self.tree
            #self.parseTreeComplete( self.tree, outname = outname )
        #f.Close()
####################################################################################################################
    def setROOTFile( self, fname ):
        self.name = None
        if os.path.exists( fname ):
            self.fname = fname
        else:
            print "no file", fname
            self.fname = None
            return None
        if '.root' in fname and not '.extra.root' in fname:
            self.fname = fname.split('.root')[0] + '.extra.root'
            if not os.path.exists( self.fname ):
                t = timer('copy file')
                shutil.copy( fname, self.fname )
                print "working copy created", self.fname
            else:
                print "working copy selected", self.fname
        return self.fname
####################################################################################################################
    def print_tree( self ):
        if not self.fname:
            print "file was not set"
            return None
        f = ROOT.TFile( self.fname )
        tree = f.Get('hitSumm')
        tree.Print()
        f.Close()
        return
####################################################################################################################
    def print_file( self ):
        if not self.fname:
            print "file was not set"
            return None
        print root_numpy.list_trees( self.fname )
        return
####################################################################################################################
    def print_branches( self, branches = [] ):
        if not self.fname:
            print "file was not set"
            return None
        self.check_branch('ID')
        for branch in branches:
            if not self.check_branch(branch):
                return None
        print ['ID'] + branches
        t = timer('access tree')
        data = root_numpy.root2array( self.fname, treename = self.DATATREE, branches = ['ID'] + branches ).view(np.recarray)
        del t
        print data, len(data)
        return
####################################################################################################################
    def check_branch( self, branch = '', outname = None ):
        if not self.fname:
            print "file was not set"
            return None
        #print root_numpy.list_branches( self.fname, treename = self.DATATREE )
        if not branch in root_numpy.list_branches( self.fname, treename = self.DATATREE ):
            return self.compute_branch( branch, outname = outname )
        else:
            return True
####################################################################################################################
    def remove_branch( self, branch_name = None, outname = None ):
        if not branch_name in root_numpy.list_branches( self.fname, treename = self.DATATREE ):
            print 'no branch', branch_name, 'in', self.fname
            return
        #return
        #if branch_name == 'ID':
        print 'removing branch', branch_name, 'from file', self.fname
        #f = ROOT.TFile(self.fname, 'update')
        #tree = f.Get( self.DATATREE )
        #tree.Print()
        
        infile = ROOT.TFile( self.fname )
        #config = infile.Get('config')
        hitSumm = infile.Get('hitSumm')
        
        hitSumm.SetBranchStatus( '*', 1 )
        hitSumm.SetBranchStatus( branch_name, 0 )
        
        if outname is None: outname = self.fname + '.tmp'
        
        outfile = ROOT.TFile( outname, 'recreate' )
        #config.CopyTree('')
        hitSumm.CopyTree('')
        outfile.Write( '', hitSumm.kOverwrite )
        outfile.Close()
        infile.Close()
    
        #print root_numpy.list_trees( self.fname )
        print root_numpy.list_trees( outname )
        print root_numpy.list_branches( outname, treename = 'hitSumm' )
        os.remove( self.fname )
        os.rename( outname, self.fname )
        return
        #else:
            #print 'only removes ID'
            #return
####################################################################################################################
    def compute_branch( self, branch = '', outname = None ):
        arr = None
        print 'computing branch', branch, 'on file', self.fname
        print root_numpy.list_branches( self.fname, treename = self.DATATREE )
        self.print_tree()
        if not outname is None:
            print 'output to', outname
        if branch in root_numpy.list_branches( self.fname, treename = self.DATATREE ):
            print 'branch is already computed; skipping'
            return True
        
        if branch == 'trackID' or branch == 'eventID':
            branches = root_numpy.list_branches( self.fname, treename = self.DATATREE )
            
            t = timer( 'compute branch ' + branch )
            data = root_numpy.root2array(self.fname, treename = self.DATATREE, branches = ['runID', 'ohdu'] ).view(np.recarray)
            runIDs = set(data['runID'])
            
            trackID = np.array( range( 1, len(data['runID']) + 1 ) )
            eventID = np.array( [ '%04d%02d%09d'%(r,o,t) for r, o, t in zip( data['runID'], data['ohdu'], trackID ) ] )
            trackID = trackID.view(np.recarray)
            eventID = eventID.view(np.recarray)
            
            if 'trackID' in branches and 'eventID' in branches:
                return True
            elif 'trackID' in branches:
                arr = np.array( zip(eventID), dtype=[('eventID', str)] ).view(np.recarray)
            elif 'eventID' in branches:
                arr = np.array( zip(trackID), dtype=[('trackID', np.int32)] ).view(np.recarray)
            else:
                arr = np.array( zip(trackID, eventID), dtype=[('trackID', np.int32), ('eventID', 'S32')] ).view(np.recarray)       

            #root_numpy.array2root( arr, self.fname, treename = self.DATATREE, mode='update' )
            #return True
        elif branch == 'aveE':
            print 'computing branch', branch, 'on file', self.fname
            t = timer('compute branch ' + branch )
            data = root_numpy.root2array(self.fname, treename = self.DATATREE, branches = ['ePix']).view(np.recarray)
            ePix = data.ePix
            n = len(ePix)
            aveE = np.array([ np.average( e[e>0] ) for e in ePix ] )
            stdE = np.array([ np.std( e[e>0] ) for e in ePix ] )
            arr = np.array( zip(aveE, stdE), dtype = [('aveE', np.float32),('stdE', np.float32)] ).view(np.recarray)
            
            #root_numpy.array2root( arr, self.fname, treename = self.DATATREE, mode='update' )
            #return True
        elif branch == 'logRedVarE':
            self.check_branch('aveE')
            t = timer('compute branch ' + branch )
            print 'computing branch', branch, 'on file', self.fname
            data = root_numpy.root2array(self.fname, treename = self.DATATREE, branches = ['aveE','stdE']).view(np.recarray)
            aveE = data['aveE']
            stdE = data['stdE']
            mask = np.logical_and( aveE>0, stdE>0 )
            print 'zeros', aveE.shape, aveE[mask].shape
            value = aveE*stdE
            value[mask] = np.log10( stdE[mask]/aveE[mask] )
            print 'zeros', value.shape, value[value==value].shape
            arr = np.array( zip(value), dtype = [(branch, np.float32)] ).view(np.recarray)
            #root_numpy.array2root( arr, self.fname, treename = self.DATATREE, mode = 'update' )
            #self.print_tree()
            #return True
        elif branch == 'T':
            print 'computing branch', branch, 'on file', self.fname
            t = timer('compute branch ' + branch )
            data = root_numpy.root2array(self.fname, treename = self.DATATREE, branches = ['xPix','yPix','ePix','nSavedPix']).view(np.recarray)
            
            T = np.zeros_like( data['nSavedPix'], dtype = np.float )
            L = np.zeros_like( data['nSavedPix'], dtype = np.float )
            stdT = np.zeros_like( data['nSavedPix'], dtype = np.float )
            stdL = np.zeros_like( data['nSavedPix'], dtype = np.float )
            rmseT = np.zeros_like( data['nSavedPix'], dtype = np.float )
            rmseL = np.zeros_like( data['nSavedPix'], dtype = np.float )
            xyBar = np.zeros_like( data['nSavedPix'], dtype = np.float )
            xBar = np.zeros_like( data['nSavedPix'], dtype = np.float )
            yBar = np.zeros_like( data['nSavedPix'], dtype = np.float )
            x2Bar = np.zeros_like( data['nSavedPix'], dtype = np.float )
            y2Bar = np.zeros_like( data['nSavedPix'], dtype = np.float )
            
            for i in range(len(data['nSavedPix'])):
                eePix = data['ePix'][i]
                mask = eePix > 0
                
                ePix = data['ePix'][i][mask]
                xPix = data['xPix'][i][mask].astype(np.float)
                yPix = data['yPix'][i][mask].astype(np.float)
                
                #print 'xPix', xPix, xPix.shape, np.max(xPix) - np.min(xPix)
                #print 'yPix', yPix, yPix.shape, np.max(yPix) - np.min(yPix)
                #xBar[i] = np.average( xPix, weights = ePix )
                #yBar[i] = np.average( yPix, weights = ePix )
                #xyBar[i] = np.average( xPix*yPix, weights = ePix )
                #x2Bar[i] = np.average( xPix**2, weights = ePix )
                #y2Bar[i] = np.average( yPix**2, weights = ePix )
                
                #vPix = np.array( zip(xPix, yPix) )
                
                if not len(ePix) > 0: continue
                C = cov( np.vstack((xPix, yPix)), ePix )
                w, v = np.linalg.eig( C )
                if w[0] < 0: 
                    w[0] *= -1
                    v[:,0] *= -1
                if w[1] < 0: 
                    w[1] *= -1
                    v[:,1] *= -1
                if w[0] > w[1]: wLong, wTrans, vLong, vTrans = w[0], w[1], v[:,0], v[:,1]
                else: wLong, wTrans, vLong, vTrans = w[1], w[0], v[:,1], v[:,0]

                lPix = np.inner( vLong, zip( xPix, yPix ) )
                #print 'lPix', lPix, lPix.shape, np.max(lPix) - np.min(lPix)
                L[i] = np.max(lPix) - np.min(lPix)
                tPix = np.inner( vTrans, zip( xPix, yPix ) )
                #print 'tPix', tPix, tPix.shape, np.max(tPix) - np.min(tPix)
                T[i] = np.max(tPix) - np.min(tPix)
                stdL[i] = std( lPix, weights = ePix )
                stdT[i] = std( tPix, weights = ePix )
                
                #print std(xPix, weights = ePix), std(yPix, weights = ePix), stdL[i], stdT[i]
                #print np.std(xPix), np.std(yPix)
                
                tMean = np.average( tPix, weights = ePix )
                lMean = np.average( lPix, weights = ePix )
                if np.max(ePix) > 0:
                    rmseT[i] = np.sqrt( np.mean( ( np.exp( -.5* (tPix - tMean)**2/stdT[i]**2) - ePix/np.max(ePix) )**2 ) ) if stdT[i] > 0 else 0
                    rmseL[i] = np.sqrt( np.mean( ( np.exp( -.5* (lPix - lMean)**2/stdL[i]**2 ) - ePix/np.max(ePix) )**2 ) ) if stdL[i] > 0 else 0
                else:
                    rmseT[i] = 0
                    rmseL[i] = 0
            
            arr = np.array( zip( T, L, stdT, stdL, rmseT, rmseL ), dtype = [
                    ('T', np.float32),
                    ('L', np.float32),
                    ('stdT', np.float32),
                    ('stdL', np.float32),
                    ('rmseT', np.float32),
                    ('rmseL', np.float32),
                ] ).view(np.recarray)
        elif branch == 'relLen':
            self.check_branch('aveE')
            t = timer('compute branch ' + branch )
            print 'computing branch', branch, 'on file', self.fname
            data = root_numpy.root2array(self.fname, treename = self.DATATREE, branches = ['T','L','stdT','stdL','rmseT','rmseL','n0','E0']).view(np.recarray)
            n0 = data['n0']
            print 'n0', data['n0'].shape
            E0 = data['E0']
            T = data['T']
            L = data['L']
            stdT = data['stdT']
            stdL = data['stdL']
            rmseT = data['rmseT']
            rmseL = data['rmseL']
            relLen = np.zeros_like( T )
            redVarT = np.zeros_like( T )
            redVarL = np.zeros_like( T )
            relVar = np.zeros_like( T )
            redn0 = np.zeros_like( T )
            logRedRmseT = np.zeros_like( T )
            relRmse = np.zeros_like( T )
            
            logEloss = np.zeros_like( T )
            
            print 'data.T', T.shape
            print 'data[T]', data['T'].shape
            relLen[ L>0 ] = T[ L>0 ]/L[ L>0 ]
            redVarT[ T>0 ] = stdT[ T>0 ]/T[ T>0 ]
            redVarL[ L>0 ] = stdL[ L>0 ]/L[ L>0 ]
            relVar[ stdL>0] = stdT[ stdL>0 ]/stdL[ stdL>0 ]
            redn0[ L*T>0] = n0[ L*T>0]/( L * T )[ L*T>0]
            logRedRmseT[ rmseT*L>0 ] = np.log10( rmseT[ rmseT*L>0 ]/L[ rmseT*L>0 ] )
            relRmse[ rmseL>0 ] = rmseT[ rmseL>0 ]/rmseL[ rmseL>0 ]
            
            logEloss[ E0>0 ] = np.log10( E0[ E0>0 ]/np.sqrt( (L[ E0>0 ]*15)**2 + 675**2 ) )
            
            print 'relVar', relVar.shape
            
            arr = np.array( zip(relLen, redVarT, redVarL, relVar, redn0, logRedRmseT, relRmse, logEloss), dtype = [
                    ('relLen', np.float32),
                    ('redVarT', np.float32),
                    ('redVarL', np.float32),
                    ('relVar', np.float32),
                    ('redn0', np.float32),
                    ('logRedRmseT', np.float32),
                    ('relRmse', np.float32),
                    ('logEloss', np.float32),
                ] ).view(np.recarray)
        else:
            print 'declaration of', branch, 'not given'
            return False
        
        if arr is None:
            print 'fatal error! arr not set', branch
        if outname is None or outname == self.fname:
            print 'updating file', self.fname
            root_numpy.array2root( arr, self.fname, treename = self.DATATREE, mode='update' )
        else:
            print 'duplicating file', self.fname
            shutil.copyfile( self.fname, outname )
            print 'updating file', outname
            root_numpy.array2root( arr, outname, treename = self.DATATREE, mode='update' )
            self.fname = outname
        
        self.print_tree()
        return True

####################################################################################################################
    def parseTree( self, tree ):
        t1 = timer("parseTree")
        print "parse tree"
        ROOT.gROOT.ProcessLine("\
            struct event_t {\
            Int_t runID;\
            Int_t ohdu;\
            Float_t E0;\
            Int_t nSat;\
            Float_t n0;\
            Float_t xVar0;\
            Float_t yVar0;\
            Float_t xBary0;\
            Float_t yBary0;\
            Int_t xPix[100000];\
            Int_t yPix[100000];\
            Float_t ePix[100000];\
            Int_t level[100000];\
            Int_t flag;\
            Int_t nSavedPix;\
            }" )
        event = ROOT.event_t()

        for var in self.vars_all: tree.SetBranchAddress( var, AddressOf( event, var ) )

        self.tracks = []
        print 'size of catalog', tree.GetEntries()
        t1.factor(text = 'entries', f = 1./tree.GetEntries() )

        #runID = []
        #ohdu = []
        #for i in xrange( tree.GetEntries() ):
            #tree.GetEntry(i)
            #runID += [event.runID]
            #ohdu += [event.ohdu]

        #print "runID", min(runID), max(runID), len(runID)

        for i in xrange( tree.GetEntries() ):
            trackID = i
            #if i > 10000: break
            #t = timer("for")
            tree.GetEntry(i)
            if event.flag == 0 and event.nSavedPix > 0:
                track = Track()
                track.id = trackID
                for var in self.vars_scalar:
                    setattr(track, var, getattr( event, var ))
                    #self.__dict__[var] += [ track.__dict__[var] ]
                for var in self.vars_list:
                    #setattr(track, var, getattr( event, var )[:event.nSavedPix] )
                    setattr(track, var, [ getattr( event, var )[pix] for pix in xrange( event.nSavedPix ) ] )
                    #self.__dict__[var] += [ track.__dict__[var] ]
                track.eigcov()
                self.tracks += [ track ]
        print 'number of tracks', len(self.tracks)
        tree.ResetBranchAddress(0)
	#tree.Delete()
        return self
####################################################################################################################
    def parseTreeComplete( self, outname = None ):
        t = timer('parseTreeComplete')
        print root_numpy.list_branches( self.fname, treename = self.DATATREE )
        
        self.check_branch('eventID')
        self.check_branch('aveE')
        self.check_branch('T')
        self.check_branch('logRedVarE', outname = outname )
        self.check_branch('relLen')
        
        #t0 = timer('root2array')
        #a = root_numpy.root2array(self.fname, treename = 'hitSumm', selection = 'nSavedPix>20 && flag==0' ).view(np.recarray)
        #del t0
        #t1 = timer('set locals')
        #localvars = self.__dict__
        #print 'types', a.dtype.names
        #for var in a.dtype.names:
            #localvars[var] = a[var]
        #del t1
        return self
####################################################################################################################
    def plotCorrelations(self):
        varsCor = ['E0', 'nSavedPix'] + self.vars_extra
        #f, ax = plt.subplots( len(varsCor), len(varsCor) )
        f = plt.figure()
        for i1, var1 in enumerate(varsCor):
            print i1, var1
            for i2, var2 in enumerate(varsCor):
                print i2, var2
                if i1 == i2: continue
                ax = f.add_subplot(len(varsCor), len(varsCor), i1+i2)
                ax.scatter( self.__dict__[var1], self.__dict__[var2] )
                ax.set_xlabel(r'%s'%var1)
                ax.set_ylabel(r'%s'%var2)
            print 'printing'
            plt.savefig('correlations.png')
            print 'done'
        plt.clf()
####################################################################################################################
    #def computeCombined(self ):
        #t = timer('computeCombined')

        #x_ = self.__dict__

        ##t0 = timer('convert array')
        ##for var in ['E0', 'T', 'L', 'stdT', 'stdL', 'rmseT', 'rmseL', 'nSavedPix', 'n0']:
            ##x_[var] = np.array( getattr(self, var) )
        ##del t0

        #t1 = timer('compute combined')
        #varsCor = ['logEloss', 'logRedRmseT', 'redn0', 'relLen', 'redVarT', 'logRedVarE' ] #, 'stdlog10E/avelog10E', 'log10(nL/nH)', 'aveE/Emax' ]

        ##x_['log10(dEdx)'] = np.log10(x_['E0']/np.sqrt(x_['L']**2*15**2 + 250**2)) #615m
        ##x_['log10(rmseT/L**2)'] = np.log10(x_['rmseT']/x_['L']**2)
        ##x_['nSavedPix*n0/(T*L)**2'] = x_['nSavedPix']*x_['n0']/(x_['L']*x_['T'] )**2
        ##x_['T/L'] = x_['T']/x_['L']
        ##x_['stdT/T'] = x_['stdT']/x_['T']
        ##x_['logRedVarE'] = 
        #del t1

        #for var in varsCor:
            #if not var in x_:
                #print var, 'not found'
                #exit(0)

        #return varsCor
####################################################################################################################
    def PCA(self, categories = None ):
        print "entering pca"
        #varsCor = self.computeCombined()

	if categories == None:
	  categories = [["muon", ['red', "$stdT < .5 and $nSavedPix > 200" ]],
		 ["muon2", ['magenta', "$stdT < .7 and $nSavedPix > 50 and $rmseT/$L < .8" ]],
		 #["single", ['blue', "nSavedPix == 9" ]],
		 ["lowE", ['yellow', "$E0 < 10**3.5 and $nSavedPix > 9" ]],
		 ["electron", ['black', "$rmseT > 50 and $stdT > .8" ]],
		 ["single", ['blue', "$stdT/$stdL > .8 and $nSavedPix > 9 and $rmseT/$L < 1.5" ]],
		 ["others", ['cyan', "$nSavedPix > 0" ]],
	      ]
	else:
	    categories = parse_tests(categories)

	print 'categories', categories
	parsed_categories = [ [key, [value[0], string_into_test2("self", value[1], list_of_vars = self.vars_scalar + self.vars_extra)]]
		      for key, value in categories ]
	print "parsed_categories", parsed_categories
	color_test = '"w"'
	for key, value in parsed_categories[::-1]:
	    color_test = '("%s" if %s else %s)'%( value[0], value[1], color_test )
	print
	print color_test
	colors2 = [ eval(color_test) for n in range(len(self.E0)) ]
        s = .6
        alpha = 1
        f = plt.figure(figsize=(20,20))
        f.subplots_adjust(hspace = 0, wspace = 0)
        x_ = self.__dict__
        axes = {}
        var_lim = {}
        index = lambda b,a: a*len(varsCor)+b+1
        print 'begin plots'
        t2 = timer('plots')
        for i1, var1 in enumerate( varsCor ):
            ax = f.add_subplot(len(varsCor), len(varsCor), index(i1,i1) )
            ax.grid(True)
            x_masked = x_[var1][ x_[var1] == x_[var1] ]
            ax.hist( x_masked, bins = int(np.sqrt(len( x_[var1] ))), log=True, histtype = 'step', lw = 2 )
            var_lim[i1] = ax.get_xlim()
            for key, value in categories:
                x_masked2 = [ x for x,color2 in zip(x_[var1], colors2) if color2 == value[0] ]
                ax.hist( x_masked2, bins = int(np.sqrt(len( x_masked ))), log=True, histtype = 'step', color = value[0] )
            ax.set_xlabel(var1)
        print 'done histograms'
        for i1, var1 in enumerate( varsCor ):
            print var1
            for i2, var2 in enumerate( varsCor[:(i1)] ):
                print 'vs', var2
                ind = index(i1,i2)
                axes[ind] = f.add_subplot(len(varsCor), len(varsCor), index(i1,i2) )
                ax = axes[ind]
                ax.grid(True,'both')
                ax.set_xlim( var_lim[i1] )
                ax.set_ylim( var_lim[i2] )
                #ax.scatter( x_[var1], x_[var2], s = s, c = colors2, alpha = alpha, linewidths = 0 )
                for key, value in categories[::-1]:
                    x_masked = [ x for x, color2 in zip(x_[var1],colors2) if color2 == value[0] ]
                    y_masked = [ x for x, color2 in zip(x_[var2],colors2) if color2 == value[0] ]
                    ax.plot( x_masked, y_masked, 'o', ms = .2, mec = value[0], mfc = value[0], label = key)#, aa = False, rasterized = True )
                ax.set_xlabel(var1)
                ax.set_ylabel(var2)
        del t2
        figname = (os.path.dirname(__file__) if os.path.dirname(__file__) != '' else '.') + '/plots/' + os.path.basename(self.fname) + '.pca_w.png'
        print os.path.dirname(__file__)
        print 'figname', figname
        t3 = timer('saving')
        plt.savefig( figname, bbox_inches='tight', pad_inches=0 )
        del t3
        print 'done'
        return figname, categories
####################################################################################################################
    def show_learn( self, user ):
        #varsCor = self.computeCombined()

        categories, done, data, data2, properties = read_categories( user, self.fname, None )
        #todo = [ False if i in done else True for i in range(len(self.nSavedPix)) ]
        filename = self.fname.split('.root')[0] + '.learn.%s.dict'%user

        #color_list = ['red', 'green', 'blue', 'black' ] 
        colormapping = {
            'm': 'red',#(1,  0,  0,  1), #red
            'n': 'green',#(.80,  .25,  .35,  1), #pink
            't': 'blue',#(.5,  .15,  .5,  1), #purple
            'e': 'pink',#(0,  0,  1,  1), #green
            'f': 'lime',#(0,  1,  1,  1), #green
            'd': 'cyan',#(0,  0,  .5, 1), #green
            's': 'olive',#(0,  1,  0,  1), #blue
            'b': 'orange',#(.75,0,  .75,1), #blue
            'l': 'brown',#(0,  .75,.75,1), #cyan
            'z': 'purple',#(.5, 0,  0,  1), #cyan
            'x': 'yellow',#(0,  0,  0,  1), #cyan
            'w': 'gray',#(0,  .5, 0,  1), #yellow
            'g': 'black',
            'p': 'purple',
            }

        s = 5.
        N = len(varsCor)

        ftrain = plt.figure(figsize=(30,30))
        flearned = plt.figure(figsize=(20,20))
        ftrain.subplots_adjust(hspace = 0, wspace = 0)
        flearned.subplots_adjust(hspace = 0, wspace = 0)
        var_lim = {}
        
        efficiency = {}
        contamination = {}
        cumFalseNegative = {}
        lastFalseNegative = {}
        cumFalsePositive = {}
        lastFalsePositive = {}
        lastN = 100
        feff = plt.figure()
        ax_eff = feff.add_subplot(211)
        ax_ctm = feff.add_subplot(212, sharex=ax_eff)
        for key, values in data.iteritems():
            cumFalseNegative[key] = []
            cumFalsePositive[key] = []
            lastFalseNegative[key] = []
            lastFalsePositive[key] = []
            falsePositive = []
            falseNegative = []
            
            for value in values:
                if 'missed' in value:
                    falsePositive += [1.]
                elif 'guess' in value:
                    falsePositive += [0.]
                    if not value['guess'] is None:
                        if value['guess'][0] == key:
                            falseNegative += [0.]
                        else:
                            falseNegative += [1.]
                        #efficiency[key] += [ [total, correct/total ] ]
                        cumFalseNegative[key] += [ [len(falseNegative), np.mean(falseNegative) ] ]
                        lastFalseNegative[key] += [ [len(falseNegative), np.mean(falseNegative[max( len(falseNegative)-lastN, 0 ):] ) ] ]
                #contamination[key] += [ [guesses, 1. - correct/guesses if guesses > 0 else 1 ] ]
                cumFalsePositive[key] += [ [len(falsePositive), np.mean(falsePositive) ] ]
                lastFalsePositive[key] += [ [len(falsePositive), np.mean(falsePositive[max( len(falsePositive)-lastN, 0 ):]  ) ] ]

            x, y = zip(*cumFalsePositive[key])
            ax_eff.semilogy( x, y, '-', color = colormapping[key], label = key )
            x, y = zip(*lastFalsePositive[key])
            ax_eff.semilogy( x, y, '.', color = colormapping[key] )
            x, y = zip(*cumFalseNegative[key])
            ax_ctm.semilogy( x, y, '-', color = colormapping[key] )
            x, y = zip(*lastFalseNegative[key])
            ax_ctm.semilogy( x, y, '.', color = colormapping[key] )
        ax_eff.grid(True)
        ax_ctm.grid(True)
        ax_eff.legend( fancybox=True, framealpha=.5)
        ax_eff.set_ylabel('false positives')
        ax_ctm.set_ylabel('false negatives')
        ax_ctm.set_xlabel('number of trained tracks')
        feff.subplots_adjust(hspace = 0)
        
        
        data = root_numpy.root2array( self.fname, treename = self.DATATREE, branches = varsCor, selection = 'n0>20&&flag==0' )
        dataCat = {}
        dataAll = {}
        for key, values in categories.iteritems():
            dataCat[key] = {}
            for trackID in values:
                track = root_numpy.root2array( self.fname, treename = self.DATATREE, branches = varsCor, start = trackID-1, stop = trackID )
                for var in varsCor:
                    try: dataCat[key][var] += [track[var][0]]
                    except: dataCat[key][var] = [track[var][0]]
                    try: dataAll[var] += [ track[var][0] ]
                    except: dataAll[var] = [ track[var][0] ]

        cov_matrix5 = {}
        rw5 = {}
        v5 = {}
        m5 = {}
        fcat = plt.figure(figsize=(N*10, 3*10))
        i = 1
        for key, values in categories.iteritems():
            ax = fcat.add_subplot( 3, N, i, sharex = ax if i>1 else None, sharey = ax if i>1 else None )
            ax.set(aspect=1)
            axNeg = fcat.add_subplot( 3, N, i+N, sharex = ax, sharey = ax )
            
            i += 1
            if len(values) > 3:
                if not key in properties: properties[key] = {}
                if not key in cov_matrix5: cov_matrix5[key] = np.zeros((N,N))
                if not key in m5: m5[key] = np.zeros(N)
                
                for i1, var1 in enumerate( varsCor ):
                    xmean = properties[key][ '<%s>'%var1 ]
                    m5[key][i1] = xmean
                    x2mean = properties[key][ '<%s.2>'%var1 ]
                    cov_matrix5[key][i1,i1] = x2mean - xmean**2
                    
                    for i2, var2 in enumerate( varsCor[:i1] ):
                        ymean = properties[key][ '<%s>'%var2 ]
                        y2mean = properties[key][ '<%s.2>'%var2 ]
                        xymean = properties[key][ '<%s.%s>'%(var1,var2) ]
                        cov_matrix5[key][i1,i2] = cov_matrix5[key][i2,i1] = xymean - xmean*ymean
                w5, v5[key] = np.linalg.eig( cov_matrix5[key] )
                rw5[key] = np.sqrt(np.abs(w5))
                
                for key2 in categories.keys():
                    if key2 == key: continue
                    datCat_new = np.inner( v5[key].T, np.array([ dataCat[key2][var] for var in varsCor ]).T )
                    ax.plot( datCat_new[0,:], datCat_new[1,:], '.', color = colormapping[key2] )
                    axNeg.plot( datCat_new[-1,:], datCat_new[-2,:], '.', color = colormapping[key2] )
                
                datCat_new = np.inner( v5[key].T, np.array([ dataCat[key][var] for var in varsCor ]).T )
                ax.plot( datCat_new[0,:], datCat_new[1,:], 'o', color = colormapping[key] )
                axNeg.plot( datCat_new[-1,:], datCat_new[-2,:], 'o', color = colormapping[key] )
                ax.set_xlabel('pc[%s]'%key)
                ax.set_ylabel('sc[%s]'%key)
                axNeg.set_xlabel('-pc[%s]'%key)
                axNeg.set_ylabel('-sc[%s]'%key)
                
        i = 1
        for key in categories.keys():
            #axMed = fcat.add_subplot( 3, N, i+2*N, sharex = ax, sharey = ax )
            i +=1
            v5Med = np.mean( [ v5[key][:,-1] ] + [ v5[key2][:,0] for key2 in categories.keys() if not key2 is key ], axis = 0 )
            for key2 in categories.keys():
                if key2 is key: continue
                datMed = np.inner( v5Med.T, np.array([ dataCat[key2][var] for var in varsCor ]).T )
                #axMed.plot( datMed, datMed, '.', color = colormapping[key2] )
            datMed = np.inner( v5Med.T, np.array([ dataCat[key][var] for var in varsCor ]).T )
            #axMed.plot( datMed, datMed, 'o', color = colormapping[key] )
            #axMed.set_xlabel('pc[%s]'%key)
            #axMed.set_ylabel('sc[%s]'%key)
            


        cat_keys = {}
        print 'plot show learn'
        t_plots = timer('plots hist')
        for i1, var1 in enumerate( varsCor ):
                ax_train = ftrain.add_subplot(len(varsCor), len(varsCor), (i1*len(varsCor))+i1+1 )
                ax_learned = flearned.add_subplot(len(varsCor), len(varsCor), (i1*len(varsCor))+i1+1 )
                x = data[var1]
                for ax in [ax_train, ax_learned]:
                    n, bins, tmp = ax.hist( x[x==x], bins = int(np.sqrt(len(x))*.5), log=True, histtype = 'step', color = 'black' )
                    ax.set_xlabel(var1)
                var_lim[var1] = ax_train.get_xlim()
                for key, values in categories.iteritems():
                    xTrain = dataCat[key][var1]
                    ax_train.hist( xTrain, bins = bins, log=True, histtype = 'step', color = colormapping[key] if key in colormapping else 'red', label = key )
                ax_train.legend(fancybox=True, framealpha=0.5)
                
                
        del t_plots
        t_plots = timer('plots scatter')
        for i1, var1 in enumerate( varsCor ):
            for i2, var2 in enumerate( varsCor[:i1] ):
                cat_keys[var1+'.'+var2] = {}

                ax_train = ftrain.add_subplot(len(varsCor), len(varsCor), (i2*len(varsCor))+i1+1 )
                ax_learned = flearned.add_subplot(len(varsCor), len(varsCor), (i2*len(varsCor))+i1+1 )
                #ax_learned.scatter( data[var1][0], data[var2][0], s = s, c = (1,1,1,0), linewidths = 0 )
                ax_learned.plot( data[var1][0], data[var2][0], 'w.' )
                for key, values in categories.iteritems():
                    ax_train.plot(dataCat[key][var1], dataCat[key][var2], '.', color = colormapping[key] )

                for ax in [ax_train, ax_learned]:
                    ax.grid(True)
                    ax.set_xlim( var_lim[var1] )
                    ax.set_ylim( var_lim[var2] )
                    ax.set_xlabel(var1)
                    ax.set_ylabel(var2)
                    ax.tick_params( labelbottom='off', labeltop='on' if i1 == 0 else 'off', labelleft='off', labelright='on' if i2 == len(varsCor)-1 else 'off' )
                    ax.xaxis.set_label_position("top")
                    ax.yaxis.set_label_position("right")

                for key, values in categories.iteritems():
                    if len(values) > 1:
                        xKey, yKey = dataCat[key][var1], dataCat[key][var2]
                        m = np.array([ np.mean(xKey), np.mean(yKey)])
                        ax_learned.text(m[0], m[1], key)
                        ax_train.text(m[0], m[1], key)
                        w, v = np.linalg.eig( np.cov(xKey,yKey) )
                        
                        cat_keys[var1+'.'+var2][key] = {'m':m, 'sigma_xy':np.array( (np.std(xKey),np.std(yKey)) )  }
                        angle = np.angle( v[0,0] + v[1,0]*1j, deg = True )
                        e = patches.Ellipse(xy = m, width = np.sqrt(abs(w[0])), height = np.sqrt(np.abs(w[1])), angle = angle, alpha = .5,
                                            facecolor = colormapping[key] if key in colormapping else 'red' )
                        ax_learned.add_artist(e)

                        angle = np.angle( v5[key][i1,0] + v5[key][i2,0]*1j, deg = True )
                        e = patches.Ellipse(xy = ( m5[key][i1], m5[key][i2] ), width = rw5[key][0], height = rw5[key][1], angle = angle, alpha = .2,
                                            facecolor = colormapping[key] if key in colormapping else 'red' )
                        ax_learned.add_artist(e)

        del t_plots

                
        #distscatter = {}
        #for scatter_key, current_scatter in cat_keys.items():
            #dist = {}
            #distscatter[scatter_key] = {}
            #for key in categories.keys():
                #if key == 'u': continue
                #dist[key] = 0
                ##dists = []
                #w = current_scatter[key]['sigma_xy']
                #n = 0
                #for key2 in categories.keys():
                    #if key2 == key or key2 == 'u': continue
                    #diff_xy = current_scatter[key]['m'] - current_scatter[key2]['m']
                    #w2 = current_scatter[key2]['sigma_xy']
                    #dist[key] += np.sqrt( diff_xy[0]**2/(w[0]*w2[0]) + diff_xy[1]**2/(w[1]*w2[1]) )
                    #n += 1
                #dist[key] /= n
                #distscatter[scatter_key][key] = dist[key]

        #for key in categories.keys():
            #if key == 'u': continue
            #print key, sorted([ [distscatter[scatter][key], scatter] for scatter in cat_keys.keys() ], key = operator.itemgetter(0))[::-1]
                    #dist[key] = np.linalg.norm(current_scatter[key]['m'] - current_scatter[key2]['m'])/( np.sqrt( current_scatter[key]['sigma']*current_scatter[key2]['sigma'] ) )
##        np.sqrt( np.inner( v[:,0], (d-mean) )**2/w[0]**2 + np.inner( v[:,1], (d-mean) )**2/w[1]**2 )


        figname = os.path.dirname(__file__) + '/plots/' + os.path.basename(filename)# + '.pca_w.png'
        t_save = timer('savefigs')
        ftrain.savefig( figname + '.train.png', bbox_inches='tight', pad_inches=0 )
        flearned.savefig( figname + '.learned.png', bbox_inches='tight', pad_inches=0 )
        feff.savefig( figname + '.efficiency.png', bbox_inches='tight', pad_inches=0 )
        fcat.savefig( figname + '.pca.png', bbox_inches='tight', pad_inches=0 )
        del t_save

        return { 'train_plot': os.path.basename(filename) + '.train.png',
                'learn_plot': os.path.basename(filename) + '.learned.png',
                'eff_plot': os.path.basename(filename) + '.efficiency.png',
                'categories': categories,
                #'todo': todo,
                #'numberOfTrained': numberOfTrained,
                }


####################################################################################################################
    def compute_gaussian( self, categories, varsCor, data = None, key = None ):
        t = timer( 'calculate dist' )
        r_ = self.__dict__

        trained = {}
        for i1, var1 in enumerate( varsCor ):
            for i2, var2 in enumerate( varsCor[:i1] ):
                for thiskey, values in categories.iteritems():
                    if not key is None and not key is thiskey: continue
                    if thiskey == 'u': continue
                    if not thiskey in trained: trained[thiskey] = {}
                    if len(values) > 3:
                        x, y = [ r_[var1][j] for j in values ], [ r_[var2][j] for j in values ]
                        w, v = np.linalg.eig( np.cov(x,y) )
                        trained[thiskey]['%s.%s'%(var1,var2)] = {
                            'eigvec': [v[:,0], v[:,1]], 
                            'eigval': [w[0], w[1]], 
                            'reigvec': [v[:,0]/w[0], v[:,1]/w[1] ],
                            'center': [ np.mean(x), np.mean(y) ]
                            }
                        #if not var1 in trained[thiskey]:
                            #trained[thiskey]['%s.%s'%(var1,var2)] = {
                                #'eigvec': [ v[:,0], v[:,1] ],
                                #'eigval': [ w[0], w[1] ],
                                #'reigvec': [ v[:,0]/w[0], v[:,1]/w[1] ],
                                #'center': [ np.mean(x) ]
                                #}                            
        del t

        return result
####################################################################################################################
    def compute_distance( self, track_id, categories, varsCor, properties = None ):
        t = timer( 'compute_dist' )

        dist_ave = dict( (k,0) for k,v in categories.iteritems() )
        dist_min = dict( (k,0) for k,v in categories.iteritems() )
        dist_5 = dict( (k,0) for k,v in categories.iteritems() )
        inv_dist = dict( (k,0) for k,v in categories.iteritems() )
        sum_inv = 0.
        N = len( varsCor )
        Norm = (N**2-N)/2

        this_track = root_numpy.root2array( self.fname, treename = self.DATATREE, branches = varsCor, start = track_id-1, stop = track_id )
        
        data = {}
        #for key, values in categories.iteritems():
            #data[key] = {}
            #for trackID in values:
                #tmp = root_numpy.root2array( self.fname, treename = self.DATATREE, branches = varsCor, start = trackID-1, stop = trackID )
                #for var in varsCor:
                    #try: data[key][var] += [ tmp[var][0] ]
                    #except: data[key][var] = [ tmp[var][0] ]
        #print 'data', data
        if properties is None: properties = {}
        
        cov_matrix5 = {}
        point5 = np.zeros(N)
        point5 = { var: this_track[var][0] for var in varsCor }
        #for i1, var1 in enumerate( varsCor ):
            #for i2, var2 in enumerate( varsCor[:i1] ):
                #for key, values in categories.iteritems():
                    #if key == 'u': continue
        for key, values in categories.iteritems():
            if key == 'u': continue
            if len(values) > 3:
                if not key in properties: properties[key] = {}
                if not key in cov_matrix5: cov_matrix5[key] = np.zeros((N,N))
                point5 = np.zeros(N)
                
                for i1, var1 in enumerate( varsCor ):
                    if '<%s>'%var1 in properties[key]:
                        xmean = properties[key][ '<%s>'%var1 ]
                    else: 
                        xmean = np.mean( data[key][var1] )

                    if '<%s.2>'%var1 in properties[key]:
                        x2mean = properties[key][ '<%s.2>'%var1 ]
                    else: 
                        x2mean = np.mean( data[key][var1]**2 )

                    #if cov_matrix5[key][i1,i1] == 0: 
                        #cov_matrix5[key][i1,i1] = x2mean - xmean**2
                        #point5[i1] = this_track[var1][0] - xmean
                        #print 'cov5', key, cov_matrix5[key]
                        #print 'point5', point5
                    cov_matrix5[key][i1,i1] = x2mean - xmean**2
                    point5[i1] = this_track[var1][0] - xmean
                    
                    for i2, var2 in enumerate( varsCor[:i1] ):
                        
                        if '<%s>'%var2 in properties[key]:
                            ymean = properties[key][ '<%s>'%var2 ]
                        else: 
                            ymean = np.mean( data[key][var2] )
                            
                        if '<%s.2>'%var2 in properties[key]:
                            y2mean = properties[key][ '<%s.2>'%var2 ]
                        else: 
                            y2mean = np.mean( data[key][var2]**2 )
                        
                        if '<%s.%s>'%(var1,var2) in properties[key]:
                            xymean = properties[key][ '<%s.%s>'%(var1,var2) ]
                        else: 
                            xymean = np.mean(  data[key][var1]*data[key][var2] )
                        
                        cov_matrix5[key][i1,i2] = cov_matrix5[key][i2,i1] = xymean - xmean*ymean
                            
                        dist2 = 0
                        covarianceMatrix = [[ x2mean-xmean**2, xymean-xmean*ymean ],[ xymean-xmean*ymean, y2mean-ymean**2 ]]
                        #print covarianceMatrix
                        w, v = np.linalg.eig( covarianceMatrix )
                        rw = np.sqrt(np.abs(w))
                        if rw[0]*rw[1] == 0 or rw[0] != rw[0] or rw[1] != rw[1]: dist2 = 0
                        else:
                            d = [ this_track[var1][0] - xmean, this_track[var2][0]-ymean ]
                            dist2 = np.sqrt(np.sum( (np.dot( d, v )/rw)**2 ))
                        #print 'dist2a', dist2, rw, xmean, v, d
                        #print ymean, np.mean(data[key][var2])
                        
                        #x, y = data[key][var1], data[key][var2]
                        #m = np.vstack((x,y))
                        #print var1, var2, key
                        #print 'cov', np.cov(x,y)
                        #print 'cov2', np.dot(m, m.T)/len(x) - np.mean(m, axis=1)[:,np.newaxis]*np.mean(m.T, axis=0)[np.newaxis,:]
                        #print 'std', np.std(x)**2
                        #xa = np.array(x)
                        #print 'std2', np.mean(xa*xa) - np.mean(xa)**2
                        
                        #dist3 = 0
                        #w, v = np.linalg.eig( np.cov(x,y) )
                        #w = np.sqrt(np.abs(w))
                        #if w[0]*w[1] == 0 or w[0] != w[0] or w[1] != w[1]:
                            #dist3 = 0
                        #else:
                            #mean = [np.mean(x), np.mean(y)]
                            #d = np.array([ this_track[var1][0], this_track[var2][0]])
                            #dist3 = np.sqrt( np.inner( v[:,0], (d-mean) )**2/w[0]**2 + np.inner( v[:,1], (d-mean) )**2/w[1]**2 )
                        #print 'dist2bb', np.sqrt( np.sum( (np.dot( d-mean, v)/w )**2 ) )
                        #print 'dist2ba', np.sqrt( np.sum( (np.dot( d-mean, v)/w ))**2 )
                        #print 'dist2b', dist3, w, np.mean(x), v, d-mean
                        #print mean[1], np.mean(y)
                        
                        dist_ave[key] += dist2/Norm
                        dist_min[key] = dist2 if dist_min[key] == 0 else min( dist_min[key], dist2 )
                        #N[key] += 1.
                        
                    #if np.min( np.abs(cov_matrix5[key]) ) > 0:
                #print 'cov5', key, cov_matrix5[key]
                #print 'point5', point5
                w5, v5 = np.linalg.eig( cov_matrix5[key] )
                #print 'w5', w5
                rw5 = np.sqrt(np.abs(w5))
                if np.min(rw5) == 0: dist_5[key] = 0
                else: 
                    dist_5[key] = np.sqrt( np.sum( (np.dot( point5, v5 )/rw5)**2 ) )
                inv_dist[key] = 1./dist_ave[key]
                sum_inv += 1./dist_ave[key]
        
        dists_ave = sorted( [ [ category, value, dist_min[category], dist_5[category], inv_dist[category]/sum_inv ] for category, value in dist_ave.iteritems() if value != 0 ], key = operator.itemgetter(1) )
        del t

        return dists_ave
####################################################################################################################
    def compute_all_distances( self, categories, varsCor ):
        t = timer( 'calculate all dists' )
        r_ = self.__dict__

        allDists = {}
        stds = {}
        N = len( varsCor )
        for key, values in categories.iteritems():
            allDists[key] = np.zeros(len(self.nSavedPix))
            stds[key] = 0
            if len(values) > 2:
                for i1, var1 in enumerate( varsCor ):
                    x = np.array( r_[var1] )
                    xKey = np.array([ r_[var1][j] for j in values ])
                    mxKey = np.mean(xKey)
                    stds[key] += np.std( [x[j] for j in categories[key] ] )

                    for i2, var2 in enumerate( varsCor[:i1] ):
                        y = np.array( r_[var2] )
                        yKey = np.array([ r_[var2][j] for j in values ])
                        myKey = np.mean(yKey)

                        m = [mxKey,myKey]
                        w_, v_ = np.linalg.eig( np.cov(xKey,yKey) )
                        w_ = np.sqrt(np.abs(w_))

                        if w_[0]*w_[1] == 0 or w_[0] != w_[0] or w_[1] != w_[1]: continue
                        w20 = 1./w_[0]**2
                        w21 = 1./w_[1]**2
                        d = np.array(zip(x,y)) - m
                        d1 = w20*np.inner( v_[:,0], d )**2
                        d2 = w21*np.inner( v_[:,1], d )**2
                        allDists[key] += np.sqrt(d1+d2)
            allDists[key] /= (N**2-N)/2
        del t
        #stdsOr = sorted( stds.items(), key = operator.itemgetter(1) )
        #print stdsOr

        #dist = dict( (k,0) for k,v in categories.iteritems() )
        #N = dict( (k,0) for k,v in categories.iteritems() )
        #for i1, var1 in enumerate( varsCor ):
            #for i2, var2 in enumerate( varsCor[:i1] ):
                #for key, values in categories.iteritems():
                    #if len(values) > 2:
                        #x, y = [ r_[var1][j] for j in values ], [ r_[var2][j] for j in values ]
                        #w, v = np.linalg.eig( np.cov(x,y) )
                        #w = np.sqrt(np.abs(w))
                        #if w[0]*w[1] == 0 or w[0] != w[0] or w[1] != w[1]:
                            #dist2 = 0
                        #else:
                            #mean = [np.mean(x), np.mean(y)]
                            #d = np.array([r_[var1][ track_id ], r_[var2][ track_id ]])
                            #dist2 = np.sqrt( np.inner( v[:,0], (d-mean) )**2/w[0]**2 + np.inner( v[:,1], (d-mean) )**2/w[1]**2 )
                        #dist[key] += dist2
                        #N[key] += 1.
        #dists = sorted( [ [ category, value/N[category] if N[category] > 0 else value ] for category, value in dist.iteritems() if value != 0 ], key = operator.itemgetter(1) )
        #del t

        return allDists
####################################################################################################################
    def plot_track_train( self, track_id, properties ):
        tnew = timer('plot_track_train')
        print 'f:plot_track_train'
        r_ = self.__dict__

        t = timer('initiate track')
        #track = Track( **dict([ [ var, r_[var][track_id] ] for var in root_numpy.list_branches( self.fname, treename = self.DATATREE ) ]) )
        data = root_numpy.root2array(self.fname, treename = self.DATATREE, start = track_id-1, stop = track_id )
        #track = Track( **dict([ [ var, r_[var][track_id] ] for var in root_numpy.list_branches( self.fname, treename = self.DATATREE ) ]) )
        track = Track( **{ var:data[var][0] for var in root_numpy.list_branches( self.fname, treename = self.DATATREE ) } )
        track.eigcov()
        #track.aveE1 = np.mean(track.ePix)
        #track.stdE1 = np.std(track.ePix)
        #track.E = sum(track.ePix)
        #print 'track aveE and stdE', track.aveE1, track.stdE1, track.aveE, track.stdE
        del t

        t = timer('plot')
        f= plt.figure(figsize=(10,10))
        axmat = plt.subplot2grid((3,3),(0,0), colspan = 2, rowspan=3)
        mat, matrix = track.plot(ax=axmat, log = True, cmap = 'rainbow', origin='lower')
        #print 'matrix=', matrix
        
        grads = np.array( [[ [ matrix[i][j+1] - matrix[i][j-1], matrix[i+1][j] - matrix[i-1][j] ] for j in range(1,len(matrix[i])-1)] for i in range(1,len(matrix)-1)] )
        grads = -grads
        
        grads2 = np.array( [[ (grads[i][j] + grads[i][j+1] + grads[i+1][j] + grads[i+1][j+1])/4 for j in range(len(grads[i])-1)] for i in range(len(grads)-1)] )
        #print 'grads', grads.shape, grads[:,:,0].shape, (range(1, grads.shape[0]+1))
        clim = mat.get_clim()
        #clim = None
        #axmat.plot( [ .5*(np.max( track.yPix )+np.min( track.yPix )) - np.min( track.yPix ) + .4*track.L*track.vLong[1], .5*(np.max( track.yPix )+np.min( track.yPix )) - np.min( track.yPix ) - .4*track.L*track.vLong[1] ], 
                #[ .5*(np.max( track.xPix )+np.min( track.xPix )) - np.min( track.xPix ) + .4*track.L*track.vLong[0], .5*(np.max( track.xPix )+np.min( track.xPix )) - np.min( track.xPix ) - .4*track.L*track.vLong[0] ]
                #, 'k-' )
        axmat.tick_params( labelbottom='on', labeltop='off', labelleft='on', labelright='off' )
        #absgrad = np.sqrt( grads[:,:,0]**2 + grads[:,:,1]**2 )
        #divgrad = np.array( [[ grads[i][j+1][0] + grads[i][j-1][0] + grads[i+1][j][1] + grads[i-1][j][1] for j in range(1,len(grads[i])-1)] for i in range(1,len(grads)-1)] )
        #rotgrad = np.array( [[ grads[i][j+1][1] - grads[i][j-1][1] - grads[i+1][j][0] + grads[i-1][j][0] for j in range(1,len(grads[i])-1)] for i in range(1,len(grads)-1)] )
        print grads.shape, len(grads), len(grads[0])
        #print range(len(grads)-1)
        #for i in range(len(grads)-1):
            #for j in range(len(grads[i])-1):
                #if np.dot(grads[i][j+1], grads[i][j]) < 0:
                    #print i, j
        a = 0
        xgrad = []
        _xgrad = []
        ygrad = []
        xygrad = []
        for i in range(len(grads)-1):
            for j in range(len(grads[i])-1):

                xdir = np.array([1,0])
                xunit = xdir
                ydir = np.array([0,1])
                yunit = ydir
                xydir = np.array([1,1])
                xyunit = xydir/np.linalg.norm(xydir)
                yxdir = np.array([1,-1])
                yxunit = xydir/np.linalg.norm(xydir)

                if np.dot( grads[i][j+1], grads[i][j] )/(np.linalg.norm(grads[i][j+1])*np.linalg.norm(grads[i][j]) ) < .5:
                    factor = -np.dot(grads[i][j+1], xdir)/np.dot(grads[i][j]-grads[i][j+1], xdir)
                    if factor < 0 or factor > 1: continue
                    _xgrad += [ ydir*(1-factor) + [ i, j ] ]
                    
                    xgrad += [ factor*grads[i][j] + (1-factor)*grads[i][j+1] ]

                if np.dot( grads[i+1][j], grads[i][j] )/(np.linalg.norm(grads[i+1][j])*np.linalg.norm(grads[i][j]) ) < .5:
                    factor = -grads[i+1][j][1]/(grads[i][j][1]-grads[i+1][j][1])
                    if factor < 0 or factor > 1: continue
                    _xgrad += [ [ i + (1-factor), j ] ]
                    
                    xgrad += [ factor*grads[i][j] + (1-factor)*grads[i+1][j] ]
                    

                if np.dot( grads[i+1][j+1], grads[i][j] )/(np.linalg.norm(grads[i+1][j+1])*np.linalg.norm(grads[i][j]) ) < .5:
                    factor = - np.dot(grads[i+1][j+1], xyunit)/np.dot( grads[i][j] - grads[i+1][j+1], xyunit )
                    if factor < 0 or factor > 1: continue
                    _xgrad += [ xydir*(1-factor) + [ i, j ] ]
                    
                    xgrad += [ factor*grads[i][j] + (1-factor)*grads[i+1][j+1] ]

                if np.dot( grads[i+1][j], grads[i][j+1] )/(np.linalg.norm(grads[i+1][j])*np.linalg.norm(grads[i][j+1]) ) < .5:
                    factor = - np.dot(grads[i+1][j], yxunit)/np.dot( grads[i][j+1] - grads[i+1][j], yxunit )
                    if factor < 0 or factor > 1: continue
                    _xgrad += [ yxdir*(1-factor) + [ i, j + 1 ] ]
                    
                    xgrad += [ factor*grads[i][j+1] + (1-factor)*grads[i+1][j] ]
                    
                #if np.dot( grads[i+1][j], grads[i][j] )/(np.linalg.norm(grads[i+1][j])*np.linalg.norm(grads[i][j]) ) < -.5:
                    
                    ##factor = .5/( np.linalg.norm(grads[i+1][j]) + np.linalg.norm(grads[i][j]) )
                    ##ygrad += [ np.array([ i + .5, j ]) + factor*(grads[i][j] + grads[i+1][j]).T ]
                    #factor = 1./( abs(grads[i+1][j][1]) + abs(grads[i][j][1]) )
                    #ygrad += [ [ i + .5 - factor*(grads[i][j][1] + grads[i+1][j][1]), j ]]

                #if np.dot( grads[i+1][j+1], grads[i][j] )/(np.linalg.norm(grads[i+1][j])*np.linalg.norm(grads[i][j+1]) ) < -.5:
                    
                    #factor = 1./( abs(grads[i+1][j][1]) + abs(grads[i][j][1]) )
                    #xygrad += [ [ i + .5 - factor*(grads[i][j][1] + grads[i+1][j+1][1]), j ]]
                    
        _xgrad = np.array(_xgrad)
        xgrad = np.array(xgrad)
        ygrad = np.array(ygrad)
        #xgrad = np.array( [
            #[ i - ( grads[i][j][1] + grads[i][j+1][1] )/(abs(grads[i][j+1][0]) + abs(grads[i][j][0])), j + .5 - .5*(grads[i][j][0] + grads[i][j+1][0])/(abs(grads[i][j+1][0]) + abs(grads[i][j][0])) ] #+ .5 + (abs(grads[i][j+1][1]) - abs(grads[i][j][1]))/(abs(grads[i][j+1][1]) + abs(grads[i][j][1])) ]
            #for i in range(len(grads)-1)            
            #for j in range(len(grads[i])-1) 
            #if np.dot( grads[i][j+1], grads[i][j] ) < -.5 ])
            #grads[i][j+1][0] >= 0 and grads[i][j][0] <= 0 and abs( grads[i][j][1] + grads[i][j+1][1] )/(abs(grads[i][j+1][0]) + abs(grads[i][j][0])) <= .5 ] )

        #xygrad = np.array( [
            #[ i +.5 - ( grads[i][j][1] + grads[i+1][j+1][1] )/(abs(grads[i+1][j+1][0]) + abs(grads[i][j][0])), j + .5 - .5*(grads[i][j][0] + grads[i+1][j+1][0])/(abs(grads[i+1][j+1][0]) + abs(grads[i][j][0])) ] #+ .5 + (abs(grads[i][j+1][1]) - abs(grads[i][j][1]))/(abs(grads[i][j+1][1]) + abs(grads[i][j][1])) ]
            #for i in range(len(grads)-1)            
            #for j in range(len(grads[i])-1) 
            #if grads[i][j+1][0] >= 0 and grads[i][j][0] <= 0 and abs( grads[i][j][1] + grads[i][j+1][1] )/(abs(grads[i][j+1][0]) + abs(grads[i][j][0])) <= .5 ] )

        #ygrad = np.array( [
            #[ i +.5 - .5*( grads[i][j][1] + grads[i+1][j][1] )/(abs(grads[i+1][j][1]) + abs(grads[i][j][1])), j - (grads[i][j][0] + grads[i+1][j][0])/(abs(grads[i+1][j][1]) + abs(grads[i][j][1])) ] #+ .5 + (abs(grads[i][j+1][1]) - abs(grads[i][j][1]))/(abs(grads[i][j+1][1]) + abs(grads[i][j][1])) ]
            #for i in range(len(grads)-1)            
            #for j in range(len(grads[i])-1) 
            #if np.dot( grads[i+1][j], grads[i][j] ) < -.5 ])
            ##if grads[i+1][j][1] >= 0 and grads[i][j][1] <= 0 and abs(grads[i][j][0] + grads[i+1][j][0])/(abs(grads[i+1][j][1]) + abs(grads[i][j][1])) <= .5 ] )
        
        #print 'xgrad', xgrad
        scale= 100
        axmat.quiver( _xgrad[:,1]+1, _xgrad[:,0]+1, xgrad[:,0], xgrad[:,1], scale=scale, color='b' )
        #axmat.scatter( ygrad[:,1]+1, ygrad[:,0]+1, c='r' )
        #print 'divgrad', divgrad
        #for i in range(len(divgrad)):
            #for j in range(len(divgrad[i])):
                #if divgrad[i][j] < -.5:
                    #c = patches.Circle( (j+2,i+2), 1./(-divgrad[i][j]), facecolor=None )
                    #axmat.add_artist(c)
        #mask = absgrad > np.mean(absgrad)#-np.std(absgrad)
        #grads[:,:][mask] = [0,0]
        #grads[:,:,1][mask] = 0
        
        #axmat.quiver( range(1, grads.shape[1]+1), range(1, grads.shape[0]+1), grads[:,:,0], grads[:,:,1], scale=scale )
        
        ax = plt.subplot2grid((3,3),(0,2))
        ax.semilogy( track.tPix - track.tMean, track.ePix, 'b.' )
        mask = np.abs(track.lPix - track.lMean)<.5
        ax.semilogy( (track.tPix - track.tMean)[mask], track.ePix[mask], 'ro' )
        y = np.linspace( -track.T/2, track.T/2, num=20 )
        ax.semilogy( y, np.exp(-.5*y**2/track.stdT**2)/(np.sqrt(2*np.pi*track.stdT))*sum(track.ePix), 'g-' )
        #ax.semilogy( y, np.ones_like(y)*track.aveE, '-' )
        ax.set_ylim( [np.min(track.ePix), np.max(track.ePix)] )
        ax.set_xlabel(r'$t$/Pixel')
        ax.set_ylabel(r'$e$/keV')

        ax = plt.subplot2grid((3,3),(1,2))
        ax.set_yscale('log')
        mystdL = std(track.lPix, weights = track.ePix)
        values = track.lPix
        weights = track.ePix
        
        mystdL2 = np.sqrt(np.average((values - np.average(values, weights = weights))**2, weights=weights))
        ax.semilogy( track.lPix - track.lMean, track.ePix, 'b.' )
        mask = np.abs(track.tPix-track.tMean)<.5
        ax.semilogy( (track.lPix - track.lMean)[mask], track.ePix[mask], 'ro' )
        x = np.linspace( -track.L/2, track.L/2, num=20 )
        ax.semilogy( x, np.exp(-.5*x**2/track.stdL**2)/(np.sqrt(2*np.pi)*track.stdL)*track.E0, 'g-' )
        ax.semilogy( x, np.exp(-.5*x**2/mystdL2**2)/(np.sqrt(2*np.pi)*mystdL2)*track.E0/15, 'm-' )
        ax.semilogy( x, np.ones_like(x)*track.aveE, 'k-' )
        ax.set_ylim( [np.min(track.ePix), np.max(track.ePix)] )
        ax.set_xlabel(r'$l$/Pixel')
        ax.set_ylabel(r'$e$/keV')


        ax = plt.subplot2grid((3,3),(2,2))
        #ax.set_xscale('log')
        #_ePix = track.ePix[:]
        #_ePix = _ePix[ _ePix > 0 ]
        #n, bins, tmp = ax.hist( np.log10( _ePix ), bins = int(np.sqrt(len(track.ePix))), log = True, histtype = 'step' )
        ##ax.semilogy( bins, log10(np.exp( -bins/track.aveE )/track.aveE*sum(track.ePix)), '-' )
        ##ax.semilogy( bins, np.log10(np.exp( -np.exp(bins)/track.stdE )/track.stdE*sum(track.ePix)), '-' )
        ##ax.semilogy( bins, -bins/track.aveE ) - np.log10(track.aveE) + sum(track.ePix), '-' )
        ##ax.semilogy( bins, -bins/track.stdE ) - np.log10(track.stdE) + sum(track.ePix), '-' )
        #ax.set_xlabel(r'log10($e$/keV)')
        #ax.set_ylabel(r'count')

        dx = 1.5
        dx = 3
        factor = .05
        x = np.arange( np.min(track.lPix), np.max(track.lPix), dx )
        y = [ np.sum( track.ePix[ np.abs( track.lPix - ix ) < dx/2 ] )/dx/1e3 for ix in x ]
        y2 = [ np.sum( track.ePix[ np.abs( track.lPix - ix ) < dx/2 ] )/dx/1e3/len(track.ePix)*10 for ix in x ]
        y2s = [ std( track.tPix[ np.logical_and( np.abs( track.lPix - ix ) < dx/2, track.ePix > factor*np.mean(track.ePix) ) ], track.ePix[ np.logical_and( np.abs( track.lPix - ix ) < dx/2, track.ePix > factor*np.mean(track.ePix) ) ] ) if np.sum(track.ePix[ np.logical_and( np.abs( track.lPix - ix ) < dx/2, track.ePix > factor*np.mean(track.ePix) ) ]) > 0 else 0 for ix in x ]
        y2mean = [ np.average( track.tPix[ np.logical_and( np.abs( track.lPix - ix ) < dx/2, track.ePix > factor*np.mean(track.ePix) ) ], weights = track.ePix[ np.logical_and( np.abs( track.lPix - ix ) < dx/2, track.ePix > factor*np.mean(track.ePix) ) ] ) if np.sum(track.ePix[ np.logical_and( np.abs( track.lPix - ix ) < dx/2, track.ePix > factor*np.mean(track.ePix) ) ]) > 0 else 0 for ix in x ]
        
        
        
        
        print 'y2mean', y2mean
        y2m = np.array( [ np.sum( track.ePix[ np.logical_and( np.abs( track.lPix - ix ) < dx/2, track.tPix <= track.tMean  )] )/dx/1e3 for ix in x ] )
        y2p = np.array( [ np.sum( track.ePix[ np.logical_and( np.abs( track.lPix - ix ) < dx/2, track.tPix >= track.tMean  )] )/dx/1e3 for ix in x ] )
        
        #xvalues = [ track.vLong[1]*( ix - track.lMean ) + track.vTrans[1]*( iy - track.tMean ) + track.yMean - np.min(track.yPix) for iy, ix in zip(y2mean, x) if abs(iy) > 0 ]
        #yvalues = [ track.vLong[0]*( ix - track.lMean ) + track.vTrans[0]*( iy - track.tMean ) + track.xMean - np.min(track.xPix) for iy, ix in zip(y2mean, x) if abs(iy) > 0 ]
        #axmat.plot( xvalues, yvalues, 'k.' )
        #for ix, iy2mean, iy2s in zip( x, y2mean, y2s ):
            #if iy2mean == 0: continue
            #xval = track.vLong[1]*( ix - track.lMean ) + track.vTrans[1]*( iy2mean - track.tMean ) + track.yMean - np.min(track.yPix)
            #yval = track.vLong[0]*( ix - track.lMean ) + track.vTrans[0]*( iy2mean - track.tMean ) + track.xMean - np.min(track.xPix)
            #axmat.plot( xval, yval, 'k.' )
            #axmat.plot( [ xval + track.vTrans[1]*iy2s/2, xval - track.vTrans[1]*iy2s/2 ], [ yval + track.vTrans[0]*iy2s/2, yval - track.vTrans[0]*iy2s/2], 'k-' )
            
        #ax.plot( x - track.tMean, y, 'r.' )
        ax.plot( x - track.tMean, y2, 'r-' )
        ax.plot( x - track.tMean, y2s, 'b-' )
        ax.plot( x - track.tMean, (y2m - y2p)**2/(y2m + y2p)**2, 'g-' )
        reldiff = (y2m - y2p)**2/(y2m + y2p)**2
        track.diffET = np.mean( reldiff[ reldiff == reldiff ] )
        track.rdiffET = np.sqrt( track.diffET )
        ax.plot( [np.min(track.lPix) - track.tMean, np.max(track.lPix) - track.tMean], [ track.rdiffET, track.rdiffET], 'g--' )
        #ax.plot( x - track.tMean, y2p, 'b-' )
        ax.plot( [np.min(track.lPix) - track.tMean, np.max(track.lPix) - track.tMean], [ 10**track.logEloss/1e3, 10**track.logEloss/1e3], 'r--' )
        ax.set_xlabel(r'$x$/Pixel')
        ax.set_ylabel(r'd$E$/d$x$(MeV/Pixel)')

        del t
        print "track_id", track_id, str(track_id)
        #base_name = os.path.basename(self.fname) + '.train' + str(track_id) + '.png'
        base_name = os.path.basename(self.fname) + '.train' + str(track_id) + '.png'
        figname = os.path.dirname(__file__) + '/plots/' + base_name
        t = timer('savefig')
        f.savefig( figname, bbox_inches='tight', pad_inches=0 )
        f.clf()
        del t
        return track, base_name, clim
####################################################################################################################
    def compute_falseMatrix(self, dists_ave, data2 ):
        actual_prob = { d[0]: (( 1 - data2['matrixFalsePositive'][dists_ave[0][0]][d[0]] ) if d[0] is dists_ave[0][0] else ( data2['matrixFalsePositive'][dists_ave[0][0]][d[0]] )) *dists_ave[0][4] / ( ( 1 - data2['matrixFalsePositive'][dists_ave[0][0]][d[0]] )*dists_ave[0][4] + sum([ data2['matrixFalsePositive'][dists_ave[0][0]][d2[0]]*d2[4] for d2 in dists_ave if not d2[0] is dists_ave[0][0] ]) ) for d in dists_ave }
        probMatrix = {}
        Norm = 0
        for d in dists_ave:
            probMatrix[ d[0] ] = {}
            for d2 in dists_ave:
                posNum = ( 1 - data2['matrixFalsePositive'][d[0]][d2[0]] ) if d2[0] is d[0] else data2['matrixFalsePositive'][d[0]][d2[0]] * data2['matrixFalseNegative'][d2[0]][d[0]]
                probMatrix[ d[0] ][ d2[0] ] = posNum * d2[4]
                Norm += probMatrix[ d[0] ][ d2[0] ]
        
        guess = [ '', 0 ]
        for d in dists_ave:
            for d2 in dists_ave:
                probMatrix[ d[0] ][ d2[0] ] /= Norm
                guess = [ d2[0], probMatrix[ d[0] ][ d2[0] ] ] if probMatrix[ d[0] ][ d2[0] ] > guess[1] else guess
        
        return probMatrix, actual_prob, guess
        
####################################################################################################################
    def train_single( self, user ):
        t = timer('train_single')
        
        categories, done, data, data2, properties = read_categories( user, self.fname, None )
        
        t0 = timer('train_single:extra')
        N = len( root_numpy.root2array(self.fname, self.DATATREE, branches = ['trackID'] ) )
        while 1:
            track_id = np.random.randint(1, N+1)
            data_ = root_numpy.root2array(self.fname, self.DATATREE, branches = ['trackID', 'flag', 'n0'], start = track_id-1, stop = track_id )
            if data_['n0'][0] > 20 and data_['flag'][0] == 0: break
        del t0
        
        dists_ave = self.compute_distance( track_id, categories, varsCor, properties )
        track, base_name, clim = self.plot_track_train( track_id, properties )

        probMatrix, actual_prob, guess = self.compute_falseMatrix( dists_ave, data2 )
        
        return {
            'figname': base_name,
            'data': data,
            'data2': data2,
            'dists': dists_ave,
            #'guess': dists_ave[0] if len(dists_ave) > 1 else None,
            'guess': guess,
            'actual_prob': actual_prob,
            'probMatrix': probMatrix,
            'categories': categories,
            'track_id': track.trackID,
            'track': track,
            'training_category': None,
            }
####################################################################################################################
####################################################################################################################
    def train_category( self, user, category ):
        categories, done, data, data2, properties = read_categories( user, self.fname, None )

        if not category in categories:
            return None
        if len(categories[category]) < 3:
            return None

        N = len( root_numpy.root2array(self.fname, self.DATATREE, branches = ['trackID'] ) )
        count = 0
        while 1:
            track_id = np.random.randint(1, N+1)
            data_ = root_numpy.root2array(self.fname, self.DATATREE, branches = ['trackID', 'flag', 'n0'], start = track_id-1, stop = track_id )
            if data_['n0'][0] > 20 and data_['flag'][0] == 0:
                count += 1
                dists_ave = self.compute_distance( track_id, categories, varsCor, properties )
                probMatrix, actual_prob, guess = self.compute_falseMatrix( dists_ave, data2 )
                if guess[0] == category: break
                if count > 100: break
                #if dists_ave[0][0] == category: break

        #track_ids = np.array(range(len( self.nSavedPix )))[todo]
        #distance_to_categories = None
        #for track_id in track_ids:
            #ave_distance_to_categories = self.compute_distance( track_id, categories, varsCor )
            #if ave_distance_to_categories[0][0] == category or ave_distance_to_categories[1][0] == category:
                #break

        track, base_name, clim = self.plot_track_train( track_id, properties )

        #probMatrix, actual_prob, guess = self.compute_falseMatrix( dists_ave, data2 )
        return {
            'figname': base_name,
            'data': data,
            'data2': data2,
            'dists': dists_ave,
            'guess': guess,
            'categories': categories,
            'actual_prob': actual_prob,            
            'probMatrix': probMatrix,
            'track_id': track.trackID,
            'track': track,
            'training_category': category,
            }
####################################################################################################################
    def reconstruct_track( self, track_id, clim, base_name ):
        r_ = self.__dict__
        mytrack = Track( **dict([ [ var, r_[var][track_id] ] for var in ["ePix","xPix","yPix"] ]) )
        E = np.sum(mytrack.ePix)
        f = plt.figure()
        ax = f.add_subplot(111)
        matrix = np.zeros( ( max(mytrack.xPix)-min(mytrack.xPix) +1, max(mytrack.yPix)-min(mytrack.yPix) +1 ) )

        results = self.gaussian_adjust_recursive( mytrack )
        print E, np.sum(mytrack.ePix)
        I = 0
        #G = lambda point: np.sum( [ res['func'](point) res in results ] )
        #diff = 0
        #for i in range(matrix.shape[0]):
            #for j in range(matrix.shape[1]):
                #diff = ( G( np.array([i + min(mytrack.xPix), j + min(mytrack.yPix)]) ) -

        for result in results:
            g = result['func']
            v = result['directions']
            m = result['center']
            w = result['std']
            I += result['integral']
            for i in range(matrix.shape[0]):
                for j in range(matrix.shape[1]):
                    matrix[i,j] += g( np.array([i + min(mytrack.xPix), j + min(mytrack.yPix)]) )
            angle = np.angle( v[1,0] + v[0,0]*1j, deg = True )
            e = patches.Ellipse(xy = (np.array(m) - [min(mytrack.xPix),min(mytrack.yPix)])[::-1], width = w[0], height = w[1], angle = angle, alpha = .1, facecolor = 'k' )
            ax.add_artist(e)
            e = patches.Ellipse(xy = (np.array(m) - [min(mytrack.xPix),min(mytrack.yPix)])[::-1], width = 2*w[0], height = 2*w[1], angle = angle, alpha = .1, facecolor = 'k' )
            ax.add_artist(e)
        print 'gaussian integral', I
        mat = ax.matshow( np.log(matrix), cmap = plt.get_cmap("rainbow") )
        mat.set_clim(clim)
        plt.colorbar( mat )
        figGauss = base_name + 'gauss.png'
        f.savefig( os.path.dirname(__file__) + '/plots/' + figGauss, bbox_inches = 'tight', pad_inches = 0 )
        f.clf()

        #results = self.gaussian_find_max( mytrack )

        #f = plt.figure()
        #ax = f.add_subplot(111)
        #matrix = np.zeros( ( max(mytrack.xPix)-min(mytrack.xPix) +1, max(mytrack.yPix)-min(mytrack.yPix) +1 ) )

        #print E, np.sum(mytrack.ePix)
        #I = 0
        #for result in zip(*results):
            #g = result[0]
            #m = result[1]
            #w = result[2]
            ##I += result[3]
            #for i in range(matrix.shape[0]):
                #for j in range(matrix.shape[1]):
                    #matrix[i,j] += g( np.array([i + min(mytrack.xPix), j + min(mytrack.yPix)]) )
            ##angle = np.angle( v[1,0] + v[0,0]*1j, deg = True )
            #e = patches.Ellipse(xy = (np.array(m) - [min(mytrack.xPix),min(mytrack.yPix)])[::-1], width = w, height = w, angle = 0, alpha = .1, facecolor = 'k' )
            #ax.add_artist(e)
            #e = patches.Ellipse(xy = (np.array(m) - [min(mytrack.xPix),min(mytrack.yPix)])[::-1], width = 2*w, height = 2*w, angle = 0, alpha = .1, facecolor = 'k' )
            #ax.add_artist(e)
        ##print matrix[matrix<0]
        #mat = ax.matshow( np.log(matrix), cmap = plt.get_cmap("rainbow") )
        #mat.set_clim(clim)
        #plt.colorbar( mat )
        #figGauss2 = base_name + 'gauss2.png'
        #f.savefig( os.path.dirname(__file__) + '/plots/' + figGauss2, bbox_inches = 'tight', pad_inches = 0 )
        #f.clf()
        #print figGauss2
        return figGauss, None
####################################################################################################################
    def gaussian_adjust_recursive( self, track, gaussians = [], prev_diff = 0, count = 0 ):
        add_parent = None
        if count == 0:
            positive_e = track.ePix > 0
            track.xPix = track.xPix[positive_e]
            track.yPix = track.yPix[positive_e]
            track.ePix = track.ePix[positive_e]
        if len(track.ePix) < 4: return add_parent

        g, center, std_, v = make_gaussian2D( track.xPix, track.yPix, track.ePix )
        center_x, center_y = center
        #center_x, center_y = np.average( track.xPix, weights = track.ePix) , np.average( track.yPix, weights = track.ePix )
        #w, v = np.linalg.eig( cov( np.vstack((track.xPix, track.yPix)), track.ePix ) )

        #std_x = std( np.inner(v[:,0], zip(track.xPix, track.yPix) ), weights = track.ePix )
        #std_y = std( np.inner(v[:,1], zip(track.xPix, track.yPix) ), weights = track.ePix )
    
        std_x, std_y = std_

        if std_x == 0 or std_y == 0: return add_parent
        if (std_x < .1 or std_y < .1) and count > 0: return add_parent

        
        if std_x > std_y:
            eigL = v[:,0]
        else:
            eigL = v[:,1]
        lPix = np.inner( eigL, zip(track.xPix, track.yPix) )
        center_l = np.average( lPix, weights = track.ePix )
        n_x, n_y = max(track.xPix)-min(track.xPix)+1, max(track.yPix)-min(track.yPix)+1

        #std_x = std_y = min(std_x, std_y)

        track_lower_half = Track()
        track_upper_half = Track()
        track_lower_half.xPix = track.xPix[lPix <= center_l ]
        track_lower_half.yPix = track.yPix[lPix <= center_l ]
        track_lower_half.ePix = track.ePix[lPix <= center_l ]
        track_upper_half.xPix = track.xPix[lPix >= center_l ]
        track_upper_half.yPix = track.yPix[lPix >= center_l ]
        track_upper_half.ePix = track.ePix[lPix >= center_l ]

        #g = lambda point: gaussian( point, np.array([center_x,center_y]), np.array([std_x,std_y]), v, integral = np.sum(track.ePix) )
        
        diff = 0
        diff_upper_half = 0
        diff_lower_half = 0
        for i in range(n_x):
            for j in range(n_y):
                x, y = i + min(track.xPix), j + min(track.yPix)
                l = np.inner(eigL, [x,y] )
                index = np.logical_and( track.xPix == x, track.yPix == y )
                index_upper_half = np.logical_and( track_upper_half.xPix == x, track_upper_half.yPix == y )
                index_lower_half = np.logical_and( track_lower_half.xPix == x, track_lower_half.yPix == y )
                G = g( np.array([x, y]) )
                E = track.ePix[ index ]
                diff += ( G - E )**2 if np.sum(index) else G**2
                diff_upper_half += ( G - E )**2 if np.sum(index_upper_half) else ( 0 if l < center_l else G**2 )
                diff_lower_half += ( G - E )**2 if np.sum(index_lower_half) else ( 0 if l > center_l else G**2 )
        #print 'single', prev_diff, diff, diff_upper_half, diff_lower_half, diff_upper_half + diff_lower_half
        #print 'relative', np.sqrt(prev_diff/sum(track.ePix)**2), np.sqrt(diff/sum(track.ePix)**2)

        #print 'single E', np.sum(track.ePix), np.sum(track_upper_half.ePix), np.sum(track_lower_half.ePix), np.sum(track_upper_half.ePix) + np.sum(track_lower_half.ePix)
        this_gaussian = {'func': g, 'std':np.array([std_x,std_y]), 'directions':v, 'center':[center_x, center_y], 'integral': np.sum(track.ePix), 'diff': diff}
        if diff_upper_half != diff_upper_half : return add_parent
        if diff_lower_half != diff_lower_half : return add_parent
        #if diff > prev_diff and count > 0 and (std_x < 1 and std_y < 1) and min(std_x,std_y)/max(std_x,std_y) > .8:
        if diff > prev_diff and count > 0 and min(std_x,std_y)/max(std_x,std_y) > .8:
            #print 'add parent'
            return add_parent
            #gaussians.append( this_gaussian )
            #return gaussians

        count += 1
        result_lower_half = self.gaussian_adjust_recursive( track_lower_half, gaussians = gaussians, prev_diff = diff_lower_half, count = count )
        #count -= 1
        if result_lower_half == add_parent:
            gaussians.append( this_gaussian )
            return gaussians
        count += 1
        result_upper_half = self.gaussian_adjust_recursive( track_upper_half, gaussians = gaussians, prev_diff = diff_upper_half, count = count )
        #count -= 1
        if result_upper_half == add_parent:
            gaussians.append( this_gaussian )
            return gaussians
        #print result_upper_half
        #gaussians.extend( result_upper_half )
        return gaussians
####################################################################################################################
    def gaussian_find_max( self, track ):
        positive_e = track.ePix > 0
        track.xPix = track.xPix[positive_e]
        track.yPix = track.yPix[positive_e]
        track.ePix = track.ePix[positive_e]

        G = []
        centers = []
        widths = []
        amps = []
        total_E = np.sum(track.ePix)
        while np.sum(track.ePix)/total_E > .15:
            index_of_max = np.argmax( track.ePix )
            x_of_max = track.xPix[index_of_max]
            y_of_max = track.yPix[index_of_max]
            radius = 2
            #min( [ x_of_max - min(track.xPix),
                            #y_of_max - min(track.yPix),
                            #max(track.xPix) - x_of_max,
                            #max(track.yPix) - y_of_max,
                            #] )
            box_around_max = (zip(*[ [x,y,e] for x,y,e in zip(track.xPix, track.yPix, track.ePix) if abs(x - x_of_max) < radius and abs(y - y_of_max) < radius ]))
            #print len(box_around_max[0])
            box_around_max[2] = np.array(box_around_max[2])
            w, v = np.linalg.eig( cov( np.vstack(( box_around_max[0], box_around_max[1])), np.array(box_around_max[2]) ) )

            std_0 = std( np.inner(v[:,0], zip(box_around_max[0], box_around_max[1]) ), weights = box_around_max[2] )
            std_1 = std( np.inner(v[:,1], zip(box_around_max[0], box_around_max[1]) ), weights = box_around_max[2] )
            std_ = min(std_0, std_1)

            e = track.ePix[index_of_max]
            E = np.sum(box_around_max[2])
            if E < 0: break
            g = lambda point: gaussian( point, np.array([x_of_max, y_of_max]), np.array([std_,std_]), integral = E)
            print g([x_of_max,y_of_max]), track.ePix[index_of_max], E
            G.append(g)
            centers.append([x_of_max, y_of_max])
            widths.append(std_)
            amps.append(track.ePix[index_of_max])

            for i in range(len(track.ePix)):
                track.ePix[i] -= g( [ track.xPix[i], track.yPix[i] ] )
            print np.sum(track.ePix)/total_E
        return G, centers, widths, amps
####################################################################################################################
    def plot_track( self, user ):
        pass
####################################################################################################################
    def learn(self, user ):
        #self.__dict__['log10(E0)'] = np.log10(np.array(self.__dict__['E0']))
        #self.__dict__['log10(dEdx)'] = np.log10(np.array(self.__dict__['E0'])/np.sqrt(np.array(self.__dict__['L'])**2*15**2 + 215**2))
        #self.__dict__['log10(rmseT/L**2)'] = np.log10(np.array(self.__dict__['rmseT'])/np.array(self.__dict__['L'])**2)
        #self.__dict__['nSavedPix*n0/(T*L)**2'] = np.array(self.__dict__['nSavedPix'])*np.array(self.__dict__['n0'])/( np.array(self.__dict__['L'])*self.__dict__['T'] )**2
        #self.__dict__['T/L'] = np.array(self.__dict__['T'])/np.array(self.__dict__['L'])
        #self.__dict__['stdT/T'] = np.array(self.__dict__['stdT'])/np.array(self.__dict__['T'])
        #varsCor = ['log10(dEdx)', 'log10(rmseT/L**2)', 'nSavedPix*n0/(T*L)**2', 'T/L', 'stdT/T' ]
        varsCor = self.computeCombined()

        #i0 = 0
        #numberOfTrained = 0
        #todo = [ True for i in range(len(self.nSavedPix)) ]
        #if os.path.exists("learn.dict"):
            #categories = eval( open("learn.dict").read() )
            #i0 = max( [ max(categories[key]) for key in categories.keys()] ) + 1
            #for key in categories.keys():
                #for value in categories[key]:
                    #todo[value] = False
                    #numberOfTrained += 1
        #else:
            #categories = {}
        #print "starting at", i0
        #print "number of trained", numberOfTrained
        categories, done, data, data2, properties = read_categories( user, self.fname, None )
        todo = [ False if i in done else True for i in range(len(self.nSavedPix)) ]


        def computeDistribution( var1, var2, key, cat ):
            dist = 0
            values = cat[key]
            r_ = self.__dict__
            x, y = [ r_[var1][j] for j in values ], [ r_[var2][j] for j in values ]
            mean = [np.mean(x), np.mean(y)]
            w, v = np.linalg.eig( np.cov(x,y) )
            w = np.sqrt(np.abs(w))
            return mean, w, v

        def computeDistance( m, w, v, i, var1, var2 ):
            d = np.array([self.__dict__[var1][i], self.__dict__[var2][i]])
            dist = np.sqrt(np.inner( v[:,0], (d-m) )**2/w[0]**2 + np.inner( v[:,1], (d-m) )**2/w[1]**2)
            return dist

        def computeFullDistance( i ):
            dist = dict( (k,0) for k,v in categories.iteritems() )
            for key, values in categories.iteritems():
                if len(values) > 2:
                    for i1, var1 in enumerate( varsCor ):
                        for i2, var2 in enumerate( varsCor[:(i1+1)] ):
                            m_, w_, v_ = computeDistribution( var1, var2, key )
                            dist2 = computeDistance( m_, w_, v_, i, var1, var2 )
                            dist[key] += dist2

        def readInput(i):
            cat = raw_input("?")
            if cat == "q":
                open("learn.dict", "w").write(str(categories))
                return False
            if cat not in categories:
                newcat = raw_input("new category:")
                if newcat in categories:
                    categories[newcat] += [i]
                else:
                    categories[newcat] = [i]
            elif cat in categories: categories[cat] += [i]
            return True

        t1 = timer('learn')
        t1.factor(text='entries', f =1./len(self.nSavedPix))
        allDists = {}
        stds = {}
        for key, values in categories.iteritems():
            allDists[key] = np.zeros(len(self.nSavedPix))
            stds[key] = 0
            if len(values) > 2:
                for i1, var1 in enumerate( varsCor ):
                    x = np.array( self.__dict__[var1] )
                    xKey = np.array([ self.__dict__[var1][j] for j in values ])
                    mxKey = np.mean(xKey)
                    stds[key] += np.std( [x[j] for j in categories[key] ] )

                    for i2, var2 in enumerate( varsCor[:(i1+1)] ):
                        if i1 == i2: continue
                        y = np.array( self.__dict__[var2] )
                        yKey = np.array([ self.__dict__[var2][j] for j in values ])
                        myKey = np.mean(yKey)

                        m = [mxKey,myKey]
                        w_, v_ = np.linalg.eig( np.cov(xKey,yKey) )
                        w_ = np.sqrt(np.abs(w_))

                        #m_, w_, v_ = computeDistribution( var1, var2, key )
                        if w_[0]*w_[1] == 0 or w_[0] != w_[0] or w_[1] != w_[1]: continue
                        w20 = 1./w_[0]**2
                        w21 = 1./w_[1]**2
                        d = np.array(zip(x,y)) - m
                        d1 = w20*np.inner( v_[:,0], d )**2
                        d2 = w21*np.inner( v_[:,1], d )**2
                        allDists[key] += d1+d2
                        #for i in range(len(self.tracks)):
                            #dist = w20*np.inner( v_[:,0], d[i] )**2 + w21*np.inner( v_[:,1], d[i] )**2
                            #if dist != dist2[i]: print "diff", dist, dist2[i]
                            #allDists[key][i] += dist
                        #d = zip(x,y) - m
                        #allDists[key] += np.sqrt( np.inner( v_[:,0], d )**2/w_[0]**2 + np.inner( v_[:,1], d )**2/w_[1]**2 )
        del(t1)
        stdsOr = sorted( stds.items(), key = operator.itemgetter(1) )
        print stdsOr

        isCat = lambda i,k: i in categories[k]
        guessCat = lambda i,k: allDists[k][i] < min([ allDists[key_][i] for key_ in allDists.keys() if key_ is not k and key_ is not 'u' ])/2 and allDists[k][i] < 10
        guess2 = lambda i,k: allDists[k][i] < min([ allDists[key_][i] for key_ in allDists.keys() if key_ is not k and key_ is not 'u' ])

        for key, values in allDists.iteritems():
            #rev = False
            print "doing key", key, stds[key]
            candidates = sorted( enumerate(list(allDists[key])), key = operator.itemgetter(1), reverse=False )
            trained = [ cand for cand in candidates if todo[cand[0]] == False and guess2(cand[0],key) and cand[0] in categories[key] ]
            guess = [ cand for cand in candidates if todo[cand[0]] == True and guess2(cand[0],key) ]
            print key, "trained"
            f, ax = plt.subplots(3,3)
            for i,cand in enumerate(trained[:9]):
                track = Track( **dict( [ [var, self.__dict__[var][cand[0]] ] for var in self.vars_listextra ] ) )
                track.plot(ax=ax[i/3,i%3], log = True, cmap = 'rainbow')
                #self.tracks[ cand[0] ].plot(ax=ax[i/3,i%3], log = True, cmap = 'rainbow')
            plt.show()
            print key, "examples"
            f, ax = plt.subplots(3,3)
            for i,cand in enumerate(guess[:9]):
                track = Track( **dict( [ [var, self.__dict__[var][cand[0]] ] for var in self.vars_listextra ] ) )
                track.plot(ax=ax[i/3,i%3], log = True, cmap = 'rainbow')
                #self.tracks[ cand[0] ].plot(ax=ax[i/3,i%3], log = True, cmap = 'rainbow')
            plt.show()
            if raw_input("train more %s: "%key) == 'y':
                for cand in candidates:
                    current = cand[0]
                    if not todo[ current ]: continue
                    currentDists = sorted( [ [k, allDists[k][current]] for k in allDists.keys() if k is not 'u' ], key = operator.itemgetter(1) )
                    if (currentDists[0][0] == key or (currentDists[0][0] == 'u' and currentDists[1][0] == key)) and currentDists[0][1] < currentDists[1][1]/3:
                        print "pretty sure", currentDists[0][0]
                    else: continue
                    f, ax = plt.subplots(1)
                    track = Track( **dict( [ [var, self.__dict__[var][current] ] for var in self.vars_listextra ] ) )
                    track.plot(ax=ax, log = True, cmap = 'rainbow')
                    plt.show()
                    if not readInput( current ): break
                    todo[ current ] = False

            #if stds[key] > 1:
                #candidates = reversed(candidates)
                #print "reversed"
                #rev = True
            for cand in candidates:
                current = cand[0]
                if not todo[ current ]: continue
                currentDists = sorted( [ [k, allDists[k][current]] for k in allDists.keys() if key is not 'u' ], key = operator.itemgetter(1) )
                #if currentDists[0][1] < currentDists[1][1]/5:
                    #print "no doubt it is a", currentDists[0][0]
                    #continue
                #if currentDists[0][1] < currentDists[1][1]/2:
                    #print "pretty sure", currentDists[0][0]
                    #continue
                #if currentDists[0][1] > currentDists[1][1]/2:
                    #print "need intervention between", currentDists[0][0], currentDists[1][0]
                    #continue
                #if currentDists[0][1] < currentDists[1][1]/2:
                    #if not rev: continue
                    #print "pretty sure"
                    #print currentDists

                if stds[key] > 1:
                    print "need training", stds[key]
                    print currentDists
                elif currentDists[0][1] > 15:
                    print "I HAVE NO IDEA"
                    print currentDists
                else: continue
                print "distance", cand[1]
                f, ax = plt.subplots(1)
                track = Track( **dict( [ [var, self.__dict__[var][current] ] for var in self.vars_listextra ] ) )
                track.plot(ax=ax, log = True, cmap = 'rainbow')
                #self.tracks[ current ].plot(ax=ax, log = True, cmap = 'rainbow')
                plt.show()
                if not readInput( current ): break
                todo[ current ] = False


        for i in range(len(self.nSavedPix)):
            if not todo[i]: continue

            print "categories of %s:" % i
            matches = {}
            #dist = {}
            dist = dict( (k,0) for k,v in categories.iteritems() )
            totalMatches = 0
            fits = dict( (k,True) for k,v in categories.iteritems() )
            for i1, var1 in enumerate( varsCor ):
                for i2, var2 in enumerate( varsCor[:(i1+1)] ):
                    for key, values in categories.iteritems():
                        upper = max( [ self.__dict__[var1][j] for j in values ] )
                        lower = min( [ self.__dict__[var1][j] for j in values ] )
                        right = max( [ self.__dict__[var2][j] for j in values ] )
                        left = min( [ self.__dict__[var2][j] for j in values ] )
                        if self.__dict__[var1][i] > lower and self.__dict__[var1][i] < upper and self.__dict__[var2][i] > left and self.__dict__[var2][i] < right:
                            if key in matches: matches[key] += 1
                            else: matches[key] = 1
                            totalMatches += 1
                        else:
                            fits[key] = False

                        if len(values) > 2:
                            #x, y = [ self.__dict__[var1][j] for j in values ], [ self.__dict__[var2][j] for j in values ]
                            #m = [np.mean(x), np.mean(y)]
                            #w, v = np.linalg.eig( np.cov(x,y) )
                            #w = np.sqrt(np.abs(w))
                            #if w[0]*w[1] == 0 or w[0] != w[0] or w[1] != w[1]:
                                #dist_ = 0
                                #dist[key] += 0
                            #else:
                                #d = np.array([self.__dict__[var1][i], self.__dict__[var2][i]])
                                #dist_ = np.sqrt(np.inner( v[:,0], (d-m) )**2/w[0]**2 + np.inner( v[:,1], (d-m) )**2/w[1]**2)
                                #dist[key] += dist_

                            m_, w_, v_ = computeDistribution( var1, var2, key )
                            if w_[0]*w_[1] == 0 or w_[0] != w_[0] or w_[1] != w_[1]: dist2 = 0
                            else: dist2 = computeDistance( m_, w_, v_, i, var1, var2 )
                            dist[key] += dist2

            tmatches = sorted( [ [k,float(v)/totalMatches] for k,v in matches.iteritems()], key = operator.itemgetter(1), reverse = True )
            dists = sorted( [ [k,v] for k,v in dist.iteritems() if v != 0 and v < 20], key = operator.itemgetter(1) )
            print 'dists', dists
            if len(dists) == 0:  print "*** ??? ***"
            elif len(dists) == 1:  print "*** %s ***"%dists[0]
            elif dists[0][1] < dists[1][1]/2: print "* %s *"%dists[0]
            #print 'matches', tmatches
            if tmatches[0][1] > 2*tmatches[1][1]: print "highly probable", tmatches[0][0]
            fullMatch = [ k for k,v in fits.iteritems() if v == True ]
            if len(fullMatch) == 1: print "***", fullMatch, "***"
            if len(fullMatch) == 0: print "*** new candidate ***"

            #for number, key in enumerate(categories.keys()): print number, key
            f, ax = plt.subplots(1)

            self.track(i).plot(ax=ax, log = True, cmap = 'rainbow')
            plt.show()
            cat = raw_input("?")
            if cat == "q":
                open("learn.dict", "w").write(str(categories))
                break
            if cat not in categories:
                newcat = raw_input("new category:")
                if newcat in categories:
                    categories[newcat] += [i]
                else:
                    categories[newcat] = [i]
            elif cat in categories: categories[cat] += [i]


        #is_muon = lambda track: track.stdT < .5 and track.nSavedPix > 200
        #maybe_muon = lambda track: track.stdT < .7 and track.nSavedPix > 50 and track.rmseT/track.L < .8 and not is_muon(track)
        #is_singlePix = lambda track: track.nSavedPix == 9
        #is_lowEnergy = lambda track: track.E0 < 10**3.5 and track.nSavedPix > 9
        #is_electron = lambda track: track.rmseT > 50 and track.stdT > 0.8
        #is_singleHit = lambda track: track.stdT/track.stdL > .8 and track.nSavedPix > 9 and track.rmseT/track.L < 1.5


        colormapping = {
            'm': 'red',#(1,  0,  0,  1), #red
            'n': 'green',#(.80,  .25,  .35,  1), #pink
            't': 'blue',#(.5,  .15,  .5,  1), #purple
            'e': 'pink',#(0,  0,  1,  1), #green
            'f': 'lime',#(0,  1,  1,  1), #green
            'd': 'cyan',#(0,  0,  .5, 1), #green
            's': 'olive',#(0,  1,  0,  1), #blue
            'b': 'orange',#(.75,0,  .75,1), #blue
            'l': 'brown',#(0,  .75,.75,1), #cyan
            'z': 'purple',#(.5, 0,  0,  1), #cyan
            'x': 'yellow',#(0,  0,  0,  1), #cyan
            'w': 'gray',#(0,  .5, 0,  1), #yellow
            'g': 'black',
            }

        s = 1.

        fguess = plt.figure(figsize=(40,40))
        ftrain = plt.figure(figsize=(40,40))
        flearned = plt.figure(figsize=(40,40))
        for i1, var1 in enumerate( varsCor ):
            for i2, var2 in enumerate( varsCor[:(i1+1)] ):

                ax_guess = fguess.add_subplot(len(varsCor), len(varsCor), (i1*len(varsCor))+i2+1 )
                ax_train = ftrain.add_subplot(len(varsCor), len(varsCor), (i1*len(varsCor))+i2+1 )
                ax_learned = flearned.add_subplot(len(varsCor), len(varsCor), (i1*len(varsCor))+i2+1 )

                print (i1*len(varsCor))+i2+1,
                print i1, var1

                x = np.array(self.__dict__[var1])
                if i1 == i2:
                    for ax in [ax_guess, ax_train, ax_learned]:
                        n, bins, tmp = ax.hist( x[x==x], bins = np.sqrt(len(x)), log=True, histtype = 'step' )
                        ax.set_xlabel(var1)


                    for key, values in categories.iteritems():
                        if key == 'u': continue
                        xKey = np.array( [ x[j] for j in range(len(self.nSavedPix)) if guess2(j,key) ] )
                        xTrain = np.array( [ x[j] for j in values ] )
                        m = np.mean(xKey)
                        ax_guess.hist( xKey, bins = bins, log=True, histtype = 'step', color = colormapping[key], label = key )
                        ax_train.hist( xTrain, bins = bins, log=True, histtype = 'step', color = colormapping[key], label = key )
                    ax_guess.legend()

                    #ax.set_ylabel(var2)
                else:
                    print '\t', i2, var2
                    y = np.array(self.__dict__[var2])
                    ax_guess.scatter( y, x, s = s, c = [ colormapping[key] for i in range(len(self.nSavedPix)) for key in categories.keys() if key != 'u' and guess2(i,key) ], linewidths = 0 )
                    ax_learned.scatter( y, x, s = s, c = (1,1,1,0), linewidths = 0 )
                    ax_train.scatter( [y[j] for key,values in categories.iteritems() for j in values if key != 'u' ], [x[j] for key,values in categories.iteritems() for j in values if key != 'u' ], s = s, c = [ colormapping[key] for key, values in categories.iteritems() for j in values if key != 'u' ], linewidths = 0 )
                    #for key, values in categories.iteritems():
                        #if key == 'u': continue
                        #print key
                        #print [y[j] for j in values ]
                        #print [x[j] for j in values ]
                        #ax_train.scatter( [y[j] for j in values ], [x[j] for j in values ], s = s, c = colormapping[key], linewidths = 0 )

                    for ax in [ax_guess, ax_train, ax_learned]:
                        ax.set_xlabel(var2)
                        ax.set_ylabel(var1)
                    for key, values in categories.iteritems():
                        #upper = max( [ self.__dict__[var1][i] for i in values ] )
                        #lower = min( [ self.__dict__[var1][i] for i in values ] )
                        #right = max( [ self.__dict__[var2][i] for i in values ] )
                        #left = min( [ self.__dict__[var2][i] for i in values ] )
                        #p = patches.Rectangle( (left, lower), right - left, upper - lower, alpha=.1, facecolor = 'r' if key is 'm' else 'm' if key is 'n' else 'k' if key is 'e' else 'b' if key is 's' else 'y' if key is 'f' else 'c' if 't' else 'w' )
                        #ax.add_patch(p)
                        if key == 'u': continue
                        if len(values) > 1:
                            xKey, yKey = np.array([ x[j] for j in values ]), np.array([ y[j] for j in values ])
                            m = [ np.mean(yKey), np.mean(xKey)]
                            ax_learned.text(m[0], m[1], key)
                            ax_guess.text(m[0], m[1], key)
                            ax_train.text(m[0], m[1], key)
                            w, v = np.linalg.eig( np.cov(yKey,xKey) )
                            angle = np.angle( v[0,0] + v[1,0]*1j, deg = True )
                            e = patches.Ellipse(xy = m, width = np.sqrt(abs(w[0])), height = np.sqrt(np.abs(w[1])), angle = angle, alpha = .5,
                                                facecolor = colormapping[key] )
                            ax_learned.add_artist(e)
        fguess.savefig('corr_guess.png')
        ftrain.savefig('corr_train.png')
        flearned.savefig('corr_learned.png')

####################################################################################################################
    def getTracks( self ):
        return self.tracks
####################################################################################################################
    def track( self, i ):
        return Track( **dict([ [var, self.__dict__[var][i] ] for var in self.vars_listextra ]) )

####################################################################################################################
####################################################################################################################



class MachineLearning:
    pass
####################################################################################################################
####################################################################################################################


class Track:
####################################################################################################################
    def __init__( self, **kwargs ):
        self.__dict__.update(kwargs)
        #for name, item in kwargs.items():
            #self.__dict__[name] = item
####################################################################################################################
    def check( self ):
        self.xPix = np.array( self.xPix )
        self.yPix = np.array( self.yPix )
        self.ePix = np.array( self.ePix )
        
        self.xPix = self.xPix[ self.ePix > 0 ].astype(np.float)
        self.yPix = self.yPix[ self.ePix > 0 ].astype(np.float)
        self.ePix = self.ePix[ self.ePix > 0 ]

        return len(self.xPix)
####################################################################################################################
    def eigcov( self ):
        self.check()
        
        #if not len(ePix) > 0: continue
        C = cov( np.vstack(( self.xPix, self.yPix)), self.ePix )
        w, v = np.linalg.eig( C )
        if w[0] < 0: 
            w[0] *= -1
            v[:,0] *= -1
        if w[1] < 0: 
            w[1] *= -1
            v[:,1] *= -1
        if w[0] > w[1]: wLong, wTrans, self.vLong, self.vTrans = w[0], w[1], v[:,0], v[:,1]
        else: wLong, wTrans, self.vLong, self.vTrans = w[1], w[0], v[:,1], v[:,0]

        self.lPix = np.inner( self.vLong, zip( self.xPix, self.yPix ) )
        self.L = np.max(self.lPix) - np.min(self.lPix)
        self.tPix = np.inner( self.vTrans, zip( self.xPix, self.yPix ) )
        self.T = np.max(self.tPix) - np.min(self.tPix)
        self.stdL = std( self.lPix, weights = self.ePix )
        self.stdT = std( self.tPix, weights = self.ePix )
        
        self.tMean = np.average( self.tPix, weights = self.ePix )
        self.lMean = np.average( self.lPix, weights = self.ePix )
        self.xMean = np.average( self.xPix, weights = self.ePix )
        self.yMean = np.average( self.yPix, weights = self.ePix )
        if np.max(self.ePix) > 0:
            self.rmseT = np.sqrt( np.mean( ( np.exp( -.5* (self.tPix - self.tMean)**2/self.stdT**2) - self.ePix/np.max(self.ePix) )**2 ) ) if self.stdT > 0 else 0
            self.rmseL = np.sqrt( np.mean( ( np.exp( -.5* (self.lPix - self.lMean)**2/self.stdL**2 ) - self.ePix/np.max(self.ePix) )**2 ) ) if self.stdL > 0 else 0
        else:
            self.rmseT = 0
            self.rmseL = 0
        ##C =  np.cov( self.xPix, self.yPix )
        #C = cov( np.vstack((self.xPix.astype('float64'), self.yPix.astype('float64'))), self.ePix )

        #w, v = np.linalg.eig( C )
        #w = np.abs(w)
        #if w[0] > w[1]: self.wLong, self.wTrans, self.vLong, self.vTrans = w[1], w[0], v[:,1], v[:,0]
        #else: self.wLong, self.wTrans, self.vLong, self.vTrans = w[0], w[1], v[:,0], v[:,1]

        #self.lPix = np.inner( self.vLong, zip( self.xPix, self.yPix ) )
        #self.L = max(self.lPix) - min(self.lPix)
        #self.lReduced = (self.lPix - min(self.lPix))/(max(self.lPix) - min(self.lPix))
        #self.tPix = np.inner( self.vTrans, zip( self.xPix, self.yPix ) )
        #self.T = max(self.tPix) - min(self.tPix)
        #self.stdL = std( self.lPix, weights = self.ePix )
        #self.stdT = std( self.tPix, weights = self.ePix )

        #tMean = np.average( self.tPix, weights = self.ePix )
        #lMean = np.average( self.lPix, weights = self.ePix )
        #if np.max(self.ePix) == 0:
            #self.rmseT = 0
            #self.rmseL = 0
        #else:
            #self.rmseT = np.sum( ( np.exp( -.5* (self.tPix - tMean)**2/self.stdT**2) - self.ePix/np.max(self.ePix) )**2 )
            #self.rmseL = np.sum( ( np.exp( -.5* (self.lPix - lMean)**2/self.stdL**2 ) - self.ePix/np.max(self.ePix) )**2 )

        ##if self.stdL < self.stdT:
            ##self.stdL, self.stdT = self.stdT, self.stdL
            ##self.L, self.T = self.T, self.L
            ##self.rmseL, self.rmseT = self.rmseT, self.rmseL

        #self.xCenter = np.average( self.xPix )
        #self.yCenter = np.average( self.yPix )

        #self.r2 = (self.xPix - self.xCenter)**2 + (self.yPix - self.yCenter)**2
        #self.phi = np.arctan2( self.yPix - self.yCenter, self.xPix - self.xCenter )
        #self.normalization = np.sum( self.r2*self.ePix )
        #self.eccentricity = {}
        #self.eccentricityRaw = {}
        #self.eccentricityRawU = {}
        #self.computeEccentricity(2)
        #self.computeEccentricity(3)
        #self.computeEccentricity(4)
        #self.ecc2 = self.eccentricity[2]
        #self.ecc3 = self.eccentricity[3]
        #self.ecc4 = self.eccentricity[4]
        #self.computeSigmasOpt2(10)
        #self.stdZ_ex10, self.stdZ0_ex10, self.stdZerr_ex10 = self.stdZ_ex, self.stdZ0_ex, self.stdZerr_ex
        #self.slope10, self.intp10, self.slopeErr10 = self.slope, self.intp, self.slopeErr
        #self.computeSigmasOpt2(30)
        #self.stdZ30, self.stdZ030, self.stdZerr30 = self.stdZ, self.stdZ0, self.stdZerr
        #self.stdZ_ex30, self.stdZ0_ex30, self.stdZerr_ex30 = self.stdZ_ex, self.stdZ0_ex, self.stdZerr_ex
        #self.slope30, self.intp30, self.slopeErr30 = self.slope, self.intp, self.slopeErr
        #self.computeSigmasOpt2(50)
        #self.stdZ50, self.stdZ050, self.stdZerr50 = self.stdZ, self.stdZ0, self.stdZerr
        #self.stdZ_ex50, self.stdZ0_ex50, self.stdZerr_ex50 = self.stdZ_ex, self.stdZ0_ex, self.stdZerr_ex
        #self.slope50, self.intp50, self.slopeErr50 = self.slope, self.intp, self.slopeErr
####################################################################################################################
    def computeEccentricity( self, n ):

        #print n,
        self.eccentricityRawU[n] = np.sum( self.r2*self.ePix * np.exp( 1.j * n * self.phi) )
        self.eccentricityRaw[n] = np.abs(self.eccentricityRawU[n])
        self.eccentricity[n] = self.eccentricityRaw[n]/self.normalization
        #self.eccentricity[n] *= np.exp( -1.j * np.angle(self.eccentricity[n]) )
        #print self.eccentricity[n],
        #self.eccentricity[n] = np.real( self.eccentricity[n] )
        #print self.eccentricity[n]
        return np.real(self.eccentricity[n])
####################################################################################################################
    #def makeBins( self, min_pixels ):
        #nbins = len(self.ePix)/min_pixels
        #bins = np.linspace( 0, 1, num=nbins)
        #if len(bins) < 2: return None
        #inds = np.digitize( self.lReduced, bins )

        #stdTBin = np.array( [ std( self.tPix[ inds==i ], weights = self.ePix[ inds==i ] ) if sum(self.ePix[ inds==i ]) > 0 else 0 for i in range(nbins) ] )
####################################################################################################################
    def computeSigmasOpt2( self, min_pixels ):
        nbins = len(self.ePix)/min_pixels
        bins = np.linspace( 0, 1, num=nbins)

        self.stdZ = 0
        self.stdZ0 = 0
        self.stdZerr = 0

        self.stdZ_ex = 0
        self.stdZ0_ex = 0
        self.stdZerr_ex = 0

        self.slope = 0
        self.intp = 0
        self.slopeErr = 0

        if len(bins) < 2: return

        dbin = bins[1] - bins[0]
        inds = np.digitize( self.lReduced, bins )

        stds = np.array( [ std( self.tPix[ inds==(i+1) ], weights = self.ePix[ inds==(i+1) ] ) if sum(self.ePix[ inds==(i+1) ]) > 0 else 0 for i in range(nbins) ] )

        z = bins + .5*dbin
        self.stdZ, self.stdZ0, r_value, p_value, self.stdZerr = stats.linregress( z, stds )
        if self.stdZ < 0:
            stds = stds[::-1]
            self.stdZ, self.stdZ0, r_value, p_value, self.stdZerr = stats.linregress( z, stds )

        diff = stds - (self.stdZ * z + self.stdZ0 )
        overSigma = diff - np.std(diff) < 0

        #if self.stdZ < 0:
            #self.stdZ0 = self.stdZ0 + self.stdZ
            #self.stdZ *= -1

        if len( stds[overSigma] ) > 2:
            self.stdZ_ex, self.stdZ0_ex, r_value, p_value, self.stdZerr_ex = stats.linregress( z[overSigma], stds[overSigma] )
            #if self.stdZ_ex < 0:
                #self.stdZ0_ex = self.stdZ0_ex + self.stdZ_ex
                #self.stdZ_ex *= -1

        bins2 = np.linspace( 0.1, 0.9, num=nbins)
        inds2 = np.digitize( self.lReduced, bins2 )
        dbin2 = bins2[1] - bins2[0]

        stds2 = np.array( [ std( self.tPix[ inds2==i ], weights = self.ePix[ inds2==i ] ) if sum(self.ePix[ inds2==i ]) > 0 else 0 for i in range(nbins) ] )
        z2 = bins2 + .5*dbin2
        self.slope, self.intp, r_value, p_value, self.slopeErr = stats.linregress( z2, stds2 )
        if self.slope < 0:
            stds2 = stds2[::-1]
            self.slope, self.intp, r_value, p_value, self.slopeErr = stats.linregress( z2, stds2 )
####################################################################################################################
    def computeEnergy( self, min_pixels ):
        nbins = len(self.ePix)/min_pixels
        bins = np.linspace( 0, 1, num=nbins)

        self.dEdx = 0
        self.dEdx0 = 0
        self.dEdxerr = 0

        self.dEdx_ex = 0
        self.dEdx0_ex = 0
        self.dEdxerr_ex = 0

        if len(bins) < 2: return

        dbin = bins[1] - bins[0]
        inds = np.digitize( self.lReduced, bins )

        stds = np.array( [ std( self.tPix[ inds==i+1 ], weights = self.ePix[ inds==i+1 ] ) if sum(self.ePix[ inds==i+1 ]) > 0 else 0 for i in range(nbins) ] )

        z = bins + .5*dbin
        self.stdZ, self.stdZ0, r_value, p_value, self.stdZerr = stats.linregress( z, stds )
        if self.stdZ < 0:
            stds = stds[::-1]
            self.stdZ, self.stdZ0, r_value, p_value, self.stdZerr = stats.linregress( z, stds )

        diff = stds - (self.stdZ * z + self.stdZ0 )
        overSigma = diff - np.std(diff) < 0

        if len( stds[overSigma] ) > 2:
            self.stdZ_ex, self.stdZ0_ex, r_value, p_value, self.stdZerr_ex = stats.linregress( z[overSigma], stds[overSigma] )
####################################################################################################################
    def draw( self, matrix = np.zeros( (2048, 4096) ), relative = True, log = False ):
        offset = (min(self.xPix), min(self.yPix)) if relative else (0, 0)
        for x,y,e in zip( self.xPix, self.yPix, self.ePix ): matrix[x - offset[0], y - offset[1]] += e if not log else (np.log(e) if e > 0 else 0)
        return matrix
####################################################################################################################
    def plot( self, ax = None, cmap = 'Greys', log = False, **kwargs ):
        ax.set_xlabel( r'$x$/Pixel' )
        ax.set_ylabel( r'$y$/Pixel' )
        #ax.pcolor( np.array(self.xPix) - min(self.xPix), np.array(self.yPix) - min(self.yPix), np.array(self.ePix), cmap = plt.get_cmap(cmap), **kwargs )
        #print 'length of xPix', len(self.xPix), len(self.yPix), len(self.ePix)
        #print self.xPix
        matrix = self.draw( np.zeros( (max(self.xPix)-min(self.xPix)+1, max(self.yPix)-min(self.yPix)+1) ), log = log )
        mat = ax.matshow( matrix, cmap = plt.get_cmap(cmap), **kwargs )
        plt.colorbar( mat )
        return mat, matrix
####################################################################################################################
    def computeSigmas( self, nbins ):
        bins = np.linspace( 0.1, 0.9, num=nbins)
        dbin = bins[1] - bins[0]
        inds = np.digitize( self.lReduced, bins )

        stdTBin = np.array( [ std( self.tPix[ inds==i ], weights = self.ePix[ inds==i ] ) if sum(self.ePix[ inds==i ]) > 0 else 0 for i in range(nbins) ] )
        slope, intercept, r_value, p_value, std_err = stats.linregress( bins + .5*dbin, stdTBin )
        if slope < 0:
            stdTBin = stdTBin[::-1]
            slope, intercept, r_value, p_value, std_err = stats.linregress( bins + .5*dbin, stdTBin )
        return stdTBin, bins, slope, intercept, std_err
####################################################################################################################
    def computeSigmasOpt( self, min_pixels ):
        nbins = len(self.ePix)/min_pixels
        bins = np.linspace( 0, 1, num=nbins)
        #bins = np.linspace( 0.1, 0.9, num=nbins)
        if len(bins) < 2:
            return  {"sigma": [0,0,0],
                "energy": [0,0,0],
                "sigmaEx": [0,0,0],
                "energyEx": [0,0,0]
                }

        dbin = bins[1] - bins[0]
        inds = np.digitize( self.lReduced, bins )

        stdTBin = np.array( [ std( self.tPix[ inds==i ], weights = self.ePix[ inds==i ] ) if sum(self.ePix[ inds==i ]) > 0 else 0 for i in range(nbins) ] )

        slope, intercept, r_value, p_value, std_err = stats.linregress( bins + .5*dbin, stdTBin )
        #if slope < 0:
            #slope = -slope
            #intercept = intercept - slope*1
            #stdTBin = stdTBin[::-1]
            #slope, intercept, r_value, p_value, std_err = stats.linregress( bins + .5*dbin, stdTBin )

        diff = stdTBin - (slope * (bins + .5*dbin) + intercept )
        overSigma = diff - np.std(diff) < 0

        stdTBinEx = None
        if len( stdTBin[overSigma] ) > 2:
            slopeEx, interceptEx, r_value, p_value, std_errEx = stats.linregress( (bins + .5*dbin)[overSigma], stdTBin[overSigma] )
            #if slopeEx < 0:
                #slopeEx, interceptEx, r_value, p_value, std_errEx = stats.linregress( (bins + .5*dbin)[overSigma], stdTBin[overSigma][::-1] )
        else:
            slopeEx, interceptEx, r_value, p_value, std_errEx = 0, 0, 0, 0, 0


        totalEnergyBin = np.array( [ np.sum( self.ePix[ inds==i ] ) for i in range(nbins) ] )

        Eslope, Eintercept, r_value, p_value, Estd_err = stats.linregress( bins + .5*dbin, totalEnergyBin )
        #if Eslope < 0:
            #totalEnergyBin = totalEnergyBin[::-1]
            #Eslope, Eintercept, r_value, p_value, Estd_err = stats.linregress( bins + .5*dbin, totalEnergyBin )

        Ediff = totalEnergyBin - (Eslope * (bins + .5*dbin) + Eintercept )
        EoverSigma = Ediff - np.std(Ediff) < 0

        totalEnergyBinEx = None
        if len( totalEnergyBin[EoverSigma] ) > 2:
            EslopeEx, EinterceptEx, r_value, p_value, Estd_errEx = stats.linregress( (bins + .5*dbin)[EoverSigma], totalEnergyBin[EoverSigma] )
            #if EslopeEx < 0:
                #totalEnergyBinEx = totalEnergyBin[::-1]
                #EslopeEx, EinterceptEx, r_value, p_value, Estd_errEx = stats.linregress( (bins + .5*dbin)[EoverSigma], totalEnergyBinEx[EoverSigma] )
        else:
            EslopeEx, EinterceptEx, r_value, p_value, Estd_errEx = 0, 0, 0, 0, 0

        return {"sigma": [slope, intercept, std_err],
                "energy": [Eslope, Eintercept, Estd_err],
                "sigmaEx": [slopeEx, interceptEx, std_errEx],
                "energyEx": [EslopeEx, EinterceptEx, Estd_errEx]}
####################################################################################################################
####################################################################################################################
####################################################################################################################

def mainHelp():
    print 'options are:'
    print '--update <catalog root file>\tupdates catalog'
    print '--test <catalog root file>\sandbox catalog'
    print '--help\t\t\t\tprints this message'

def mainUpdate( fname ):
    print 'update for catalog file', fname
    if 'runID' in fname:
        outname = fname.split('.root')[0] + 'Extra.root'
        if os.path.exists( outname ):
            print 'file', outname, 'already done; not overwriting'
            exit(0)
        print 'update for catalog file', fname, 'for', outname
        catalog = Catalog( fname )
        catalog.parseTreeComplete( outname = outname )
        return True
    else:
        print 'catalog seems incorrect'
        return False

if __name__ == "__main__":
    if '--update' in sys.argv:
        print "update"
        mainUpdate( sys.argv[-1] )
        exit(0)
    elif '--split' in sys.argv:
        print "split"
        catalog = Catalog2()
        #catalog.listRunID( sys.argv[-1], nEntries = 1 )
        catalog.splitCatalog( sys.argv[-1] )
        exit(0)
    elif '--addid' in sys.argv:
        print "addid"
        catalog = Catalog2()
        #catalog.listRunID( sys.argv[-1], nEntries = 1 )
        catalog.addID( sys.argv[-1] )
        exit(0)
    elif '--analyze' in sys.argv:
        print "analyze"
        catalog = Catalog(sys.argv[-1], complete = True)
        catalog.PCA()
        exit(0)
    elif '--learn' in sys.argv:
        print "learn"
        catalog = Catalog(sys.argv[-1], complete = True)
        catalog.learn()
        exit(0)
    else:
        print 'no recognized option given'
        mainHelp()
        exit(0)
    exit(0)
