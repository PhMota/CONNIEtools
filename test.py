#!/usr/bin/env python
import ROOT
from ROOT import AddressOf
#import ctypes
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
import matplotlib.colors
import matplotlib.patches as patches
import numpy as np
#import array
#import copy
from scipy import stats
import time, sys, os
import operator

#plt.rc('text', usetex=True)

iaddattr = lambda object, name, value: setattr( object, name, getattr( object, name ) + value )


#class timer{
    #private:
        #string name;
        #double start;
    #public:
        #void ~timer( name_ ){
            #name = name_;
            #start = ...;
            #}
        
#}

class timer:
    def __init__(self, name = ''):
        self.name = name
        self.start = time.clock()
        self.text = ''
    def factor( self, text = '', f = None ):
        self.text = text
        self.f = f
    def __del__(self):
        elapsed = time.clock() - self.start
        print 'timer[%s]: %.2gs'%(self.name, elapsed)
        if not self.text == '': print '\t%s: %.3gs'%(self.text, self.f*elapsed)

def std( x, weights ):
    weights[weights!=weights] = 0
    std2 = np.average( x**2, weights = weights ) - np.average( x, weights = weights )**2
    return np.sqrt( std2 ) if std2 > 0 else np.sqrt( np.abs(std2) )

def cov(m, a):
    ddof = len(a)-1
    v1 = np.sum(a) # N
    v2 = np.sum(a**2) # N
    m -= np.sum(m * a, axis=1, keepdims=True) / v1 # m = m - (sum m)/N
    return np.dot(m * a, m.T) * v1 / (v1**2 - ddof * v2) # m*m.T * N/(N**2 - (N-1)*N)

class Catalog2:
    def listRunID( self, fname ):
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
        
        #print runID.items()
        print "runID", min(runID.keys()), max(runID.keys()), len(runID.keys())
        
        nEntries = 100
        step = 0
        minIndex = 0
        while minIndex < max(runID.keys()):
            cutfname = fname.split('.root')[0] + 'cut%03i.root' % step
            cuttfile = ROOT.TFile.Open( cutfname, "RECREATE" )
            config.CopyTree("")
            minIndex = min(runID) + step*nEntries
            maxIndex = min(runID) + (step+1)*nEntries
            tree.CopyTree("flag == 0 && nSavedPix > 0 && runID >= %i && runID < %i" % ( minIndex, maxIndex ) )
            cuttfile.Write()
            cuttfile.Close()
            step += 1
            minIndex = min(runID) + step*nEntries
        tfile.Close()
        print step
        return step
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
    def __init__( self, fname = None, opt = False, complete = False ):
        self.tree = None
        self.config = ROOT.TTree
        #if fname == None:
            #self.fname = '/home/mota/Projects/Connie/damic/damicViewer/catalog_55Fe.root'
        #else:
        self.fname = fname
        #tree = None
        self.listvars = ['xPix','yPix','ePix','level']
        #self.scalarvars = ['runID','flag','nSavedPix','ohdu','E0','n0','nSat','xVar0','yVar0','xBary0','yBary0']
        self.scalarvars = ['flag', 'runID','nSavedPix','ohdu','E0','n0']
        self.listExtra = [
            'L', 'T',
            'stdL', 'stdT',
            'errL', 'errT',
            'ecc2', 'ecc3', 'ecc4',
            'slope30', 'intp30', 'slopeErr30',
            'slope50', 'intp50', 'slopeErr50',
            'stdZ10', 'stdZ010', 'stdZerr10',
            'stdZ_ex10', 'stdZ0_ex10', 'stdZerr_ex10',
            'stdZ30', 'stdZ030', 'stdZerr30',
            'stdZ_ex30', 'stdZ0_ex30', 'stdZerr_ex30',
            #'slope30', 'intp30', 'err30',
            #'Eslope30', 'Eintp30', 'Eerr30',
            #'slopeEx30', 'intpEx30', 'errEx30',
            #'EslopeEx30', 'EintpEx30', 'EerrEx30',
            'stdZ50', 'stdZ050', 'stdZerr50',
            'stdZ_ex50', 'stdZ0_ex50', 'stdZerr_ex50',
            #'slope50', 'intp50', 'err50',
            #'Eslope50', 'Eintp50', 'Eerr50',
            #'slopeEx50', 'intpEx50', 'errEx50',
            #'EslopeEx50', 'EintpEx50', 'EerrEx50',
                ]

        self.vars = self.scalarvars + self.listvars
        
        #for var in self.vars: self.__dict__[var] = []

        #if not fname == None:
        print fname
        f = ROOT.TFile( fname )
        #self.config = f.Get("config").CloneTree()
        #self.config = f.Get("config").CloneTree()
        #print self.config
        self.tree = f.Get("hitSumm")
        if not complete:
            self.parseTree( self.tree )
        else:
            print self.tree
            self.parseTreeComplete( self.tree )
        f.Close()
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
        
        for var in self.vars: tree.SetBranchAddress( var, AddressOf( event, var ) )
            
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
                for var in self.scalarvars:
                    setattr(track, var, getattr( event, var ))
                    #self.__dict__[var] += [ track.__dict__[var] ]
                for var in self.listvars:
                    setattr(track, var, [ getattr( event, var )[pix] for pix in xrange( event.nSavedPix ) ] )
                    #self.__dict__[var] += [ track.__dict__[var] ]
                track.eigcov()
                self.tracks += [ track ]
        print 'number of tracks', len(self.tracks)
        return self
####################################################################################################################        
    def write( self ):
        f = ROOT.TFile(self.fname, 'update')
        tree = f.Get("hitSumm")
        print tree.GetEntries()
        name = self.fname.split('.root')[0] + 'Extra.root'
        print 'creating catalog', name
        
        output = ROOT.TFile.Open(name, 'recreate')
        copyTree = tree.CopyTree("flag == 0 && nSavedPix > 0")
        copyTree.SetName("hitSumm")
        
        if not copyTree.GetEntries() == len( self.tracks ):
            print 'different sizes', copyTree.GetEntries(), len( self.tracks )
            exit(0)
        
        ROOT.gROOT.ProcessLine("\
            struct event2_t {\
            Float_t T;\
            Float_t L;\
            Float_t stdT;\
            Float_t stdL;\
            Float_t errT;\
            Float_t errL;\
            Float_t ecc2;\
            Float_t ecc3;\
            Float_t ecc4;\
            Float_t slope30;\
            Float_t intp30;\
            Float_t slopeErr30;\
            Float_t slope50;\
            Float_t intp50;\
            Float_t slopeErr50;\
            Float_t stdZ10;\
            Float_t stdZ010;\
            Float_t stdZerr10;\
            Float_t stdZ_ex10;\
            Float_t stdZ0_ex10;\
            Float_t stdZerr_ex10;\
            Float_t stdZ30;\
            Float_t stdZ030;\
            Float_t stdZerr30;\
            Float_t stdZ50;\
            Float_t stdZ050;\
            Float_t stdZerr50;\
            Float_t stdZ_ex30;\
            Float_t stdZ0_ex30;\
            Float_t stdZerr_ex30;\
            Float_t stdZ_ex50;\
            Float_t stdZ0_ex50;\
            Float_t stdZerr_ex50;\
            }" )
        event = ROOT.event2_t()
        
        brch = {}
        for var in self.listExtra: brch[var] = copyTree.Branch(var, AddressOf( event, var ), '%s/F'%var )
        for track in self.tracks:
            for var in self.listExtra: setattr(event, var, getattr(track,var) )
            for var in self.listExtra: brch[var].Fill()
        
        f.Close()
        output.Write()
        output.Close()
####################################################################################################################
    def parseTreeComplete( self, tree ):
        t1 = timer("parseTree")
        print "parse tree"
        ROOT.gROOT.ProcessLine("\
            struct event_t {\
            Int_t xPix[100000];\
            Int_t yPix[100000];\
            Int_t level[100000];\
            Float_t ePix[100000];\
            Int_t flag;\
            Int_t runID;\
            Int_t ohdu;\
            Float_t E0;\
            Float_t n0;\
            Int_t nSavedPix;\
            Float_t T;\
            Float_t L;\
            Float_t stdT;\
            Float_t stdL;\
            Float_t errT;\
            Float_t errL;\
            Float_t ecc2;\
            Float_t ecc3;\
            Float_t ecc4;\
            Float_t slope30;\
            Float_t intp30;\
            Float_t slopeErr30;\
            Float_t slope50;\
            Float_t intp50;\
            Float_t slopeErr50;\
            Float_t stdZ10;\
            Float_t stdZ010;\
            Float_t stdZerr10;\
            Float_t stdZ_ex10;\
            Float_t stdZ0_ex10;\
            Float_t stdZerr_ex10;\
            Float_t stdZ30;\
            Float_t stdZ030;\
            Float_t stdZerr30;\
            Float_t stdZ50;\
            Float_t stdZ050;\
            Float_t stdZerr50;\
            Float_t stdZ_ex30;\
            Float_t stdZ0_ex30;\
            Float_t stdZerr_ex30;\
            Float_t stdZ_ex50;\
            Float_t stdZ0_ex50;\
            Float_t stdZerr_ex50;\
            }" )
        event = ROOT.event_t()
        self.varsExtra = self.scalarvars + self.listExtra + self.listvars
        for var in self.varsExtra: 
            tree.SetBranchAddress( var, AddressOf( event, var ) )
            self.__dict__[var] = []
            
        self.tracks = []
        t1.factor(text = 'entries', f = 1./tree.GetEntries() )
        
        for i in xrange( tree.GetEntries() ):
            trackID = i
            #if i > 10000: break
            #t = timer("for")
            tree.GetEntry(i)
            if event.nSavedPix > 20: # and event.ohdu == 6:
                track = Track()
                track.id = trackID
                for var in self.varsExtra:
                    setattr(track, var, getattr( event, var ) )
                    self.__dict__[var] += [ getattr( event, var ) ]
                for var in self.listvars:
                    setattr(track, var, [ getattr( event, var )[pix] for pix in xrange( event.nSavedPix ) ] )
                self.tracks += [ track ]

        return self
####################################################################################################################
    def plotCorrelations(self):
        varsCor = ['E0', 'nSavedPix'] + self.listExtra
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
    def PCA(self):
        #self.__dict__['log10(E0)'] = np.log10(np.array(self.__dict__['E0']))
        self.__dict__['log10(dEdx)'] = np.log10(np.array(self.__dict__['E0'])/np.sqrt(np.array(self.__dict__['L'])**2*15**2 + 215**2))
        self.__dict__['log10(errT/L**2)'] = np.log10(np.array(self.__dict__['errT'])/np.array(self.__dict__['L'])**2)
        self.__dict__['nSavedPix*n0/(T*L)**2'] = np.array(self.__dict__['nSavedPix'])*np.array(self.__dict__['n0'])/( np.array(self.__dict__['L'])*self.__dict__['T'] )**2
        self.__dict__['T/L'] = (np.array(self.__dict__['T'])/np.array(self.__dict__['L']))
        self.__dict__['stdT/T'] = np.array(self.__dict__['stdT'])/np.array(self.__dict__['T'])
        varsCor = ['log10(dEdx)', 'log10(errT/L**2)', 'nSavedPix*n0/(T*L)**2', 'T/L', 'stdT/T' ]# + [var for var in self.listExtra if '10' not in var and '30' not in var and '50' not in var]
        M = np.array([ self.__dict__[var] # if (var != 'dEdx' and var != 'nSavedPix' and var != 'L' and var != 'T' and var != 'stdL' and var != 'stdT' ) 
                                            # ) else np.log10(self.__dict__[var]) 
                                              for var in varsCor ])
        print M.shape
        #cov = np.cov(M)
        #print cov.shape
        ##print cov
        #cov[cov!=cov] = 0
        #w, v = np.linalg.eig( cov )
        #print v.shape
        ##pc0 = np.log10(np.dot( v[:,0], M ))
        #pc0 = np.dot( v[:,0], M )
        #pc1 = np.dot( v[:,1], M )
        ##pc1 = np.log10( pc1 - np.min(pc1))
        #pc2 = np.dot( v[:,2], M )
        ##pc2 = np.log10( pc2 - np.min(pc2))
        #pc3 = np.dot( v[:,3], M )
        ##pc3 = np.log10( pc0 - np.min(pc3))
        #print 'eigenvalues', w
        #print zip( varsCor, v[:,0] )
        #print zip( varsCor, v[:,1] )
        #print v[:,0]
        #print v[:,1]
        #print v[:,2]
        #print pc0.shape
        is_muon = lambda track: track.stdT < .5 and track.nSavedPix > 200
        maybe_muon = lambda track: track.stdT < .7 and track.nSavedPix > 50 and track.errT/track.L < .8 and not is_muon(track)
        is_singlePix = lambda track: track.nSavedPix == 9
        is_lowEnergy = lambda track: track.E0 < 10**3.5 and track.nSavedPix > 9
        is_electron = lambda track: track.errT > 50 and track.stdT > 0.8
        is_singleHit = lambda track: track.stdT/track.stdL > .8 and track.nSavedPix > 9 and track.errT/track.L < 1.5
        colors = [ 
            "r" if is_muon(track) else 
            ("m" if maybe_muon(track) else
            ("y" if is_lowEnergy(track) else 
            ("k" if is_electron(track) else 
            ( "b" if is_singleHit(track) else 
             "c"))))
            for track in self.tracks ]
        s = .6
        alpha = 1
        f = plt.figure(figsize=(25,25))
        for i1, var1 in enumerate( varsCor ):
            for i2, var2 in enumerate( varsCor[:(i1+1)] ):
                ax = f.add_subplot(len(varsCor), len(varsCor), (i1*len(varsCor))+i2+1 )
                print (i1*len(varsCor))+i2+1,
                print i1, var1
                if i1 == i2:
                    ax.hist( np.array(self.__dict__[var1])[np.array(self.__dict__[var1])==np.array(self.__dict__[var1])], bins = np.sqrt(len(self.__dict__[var1])), log=True )
                    ax.set_xlabel(var1)
                    #ax.set_ylabel(var2)
                else:
                    print '\t', i2, var2
                    ax.scatter( self.__dict__[var2], self.__dict__[var1], s = s, c = colors, alpha = alpha, linewidths = 0 )
                    ax.set_xlabel(var2)
                    ax.set_ylabel(var1)
        plt.savefig('corr_w0.png')
        #ax = f.add_subplot(331)
        #ax.scatter( self.__dict__['log10(dEdx)'], self.__dict__['log10(errT/L**2)'], s = s, c = colors, alpha = alpha, linewidths = 0 )
        #ax.set_xlabel('log10(dEdx)')
        #ax.set_ylabel('log10(errT/L**2)')
        ##ax.set_ylim((-2*np.std(pc0)+np.mean(pc0),2*np.std(pc0)+np.mean(pc0)))
        ##ax.set_xlim((-2*np.std(pc1)+np.mean(pc1),2*np.std(pc1)+np.mean(pc1)))
        
        #ax = f.add_subplot(332)
        #ax.scatter( pc2, pc0, s=s, c = colors, alpha = alpha, linewidths = 0)
        ##ax.set_ylim( (-np.std(pc0)+np.mean(pc0),np.std(pc0)+np.mean(pc0)) )
        ##ax.set_xlim( (-np.std(pc2)+np.mean(pc2),np.std(pc2)+np.mean(pc2)) )
        
        #ax = f.add_subplot(333)
        #ax.scatter( pc3, pc0, s=s, c = colors, alpha = alpha, linewidths = 0 )
        
        #ax = f.add_subplot(335)
        #ax.scatter( pc2, pc1, s=s, c = colors, alpha = alpha, linewidths = 0 )
        #ax = f.add_subplot(336)
        #ax.scatter( pc3, pc1, s=s, c = colors, alpha = alpha, linewidths = 0 )
        
        #ax = f.add_subplot(337)
        ##ax.title( "red: muons, yellow: single pixel, green: electrons, blue: circular, magenta: unclassified" )
        
        #ax = f.add_subplot(339)
        #ax.scatter( pc3, pc2, s=s, c = colors, alpha = alpha, linewidths = 0 )
        #plt.savefig('pca.png')
        #plt.show()
        #for i1, var1 in enumerate(varsCor):
            #print i1, var1
            #for i2, var2 in enumerate(varsCor):
                #print i2, var2
                #if i1 == i2: continue
                #ax = f.add_subplot(len(varsCor), len(varsCor), i1+i2)
                #ax.scatter( self.__dict__[var1], self.__dict__[var2] )
                #ax.set_xlabel(r'%s'%var1)
                #ax.set_ylabel(r'%s'%var2)
            #print 'printing'
            #plt.savefig('correlations.png')
            #print 'done'
        #plt.clf()
####################################################################################################################
    def learn(self):
        self.__dict__['log10(E0)'] = np.log10(np.array(self.__dict__['E0']))
        self.__dict__['log10(dEdx)'] = np.log10(np.array(self.__dict__['E0'])/np.sqrt(np.array(self.__dict__['L'])**2*15**2 + 215**2))
        self.__dict__['log10(errT/L**2)'] = np.log10(np.array(self.__dict__['errT'])/np.array(self.__dict__['L'])**2)
        self.__dict__['nSavedPix*n0/(T*L)**2'] = np.array(self.__dict__['nSavedPix'])*np.array(self.__dict__['n0'])/( np.array(self.__dict__['L'])*self.__dict__['T'] )**2
        self.__dict__['T/L'] = np.array(self.__dict__['T'])/np.array(self.__dict__['L'])
        self.__dict__['stdT/T'] = np.array(self.__dict__['stdT'])/np.array(self.__dict__['T'])
        varsCor = ['log10(dEdx)', 'log10(errT/L**2)', 'nSavedPix*n0/(T*L)**2', 'T/L', 'stdT/T' ]
        
        i0 = 0
        numberOfTrained = 0
        todo = [ True for i in range(len(self.tracks)) ] 
        if os.path.exists("learn.dict"):
            categories = eval( open("learn.dict").read() )
            i0 = max( [ max(categories[key]) for key in categories.keys()] ) + 1
            for key in categories.keys():
                for value in categories[key]: 
                    todo[value] = False
                    numberOfTrained += 1
        else:
            categories = {}
        print "starting at", i0
        print "number of trained", numberOfTrained
        
        def computeDistribution( var1, var2, key ):
            dist = 0
            values = categories[key]
            x, y = [ self.__dict__[var1][j] for j in values ], [ self.__dict__[var2][j] for j in values ]
            m = [np.mean(x), np.mean(y)]
            w, v = np.linalg.eig( np.cov(x,y) )
            w = np.sqrt(np.abs(w))
            return m, w, v
        
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
        t1.factor(text='entries', f =1./len(self.tracks))
        allDists = {}
        stds = {}
        for key, values in categories.iteritems():
            allDists[key] = np.zeros(len(self.tracks))
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
                self.tracks[ cand[0] ].plot(ax=ax[i/3,i%3], log = True, cmap = 'rainbow')
            plt.show()
            print key, "examples"
            f, ax = plt.subplots(3,3)
            for i,cand in enumerate(guess[:9]):
                self.tracks[ cand[0] ].plot(ax=ax[i/3,i%3], log = True, cmap = 'rainbow')
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
                    self.tracks[ current ].plot(ax=ax, log = True, cmap = 'rainbow')
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
                self.tracks[ current ].plot(ax=ax, log = True, cmap = 'rainbow')
                plt.show()
                if not readInput( current ): break
                todo[ current ] = False

        
        for i in range(len(self.tracks)):
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
            self.tracks[i].plot(ax=ax, log = True, cmap = 'rainbow')
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
        #maybe_muon = lambda track: track.stdT < .7 and track.nSavedPix > 50 and track.errT/track.L < .8 and not is_muon(track)
        #is_singlePix = lambda track: track.nSavedPix == 9
        #is_lowEnergy = lambda track: track.E0 < 10**3.5 and track.nSavedPix > 9
        #is_electron = lambda track: track.errT > 50 and track.stdT > 0.8
        #is_singleHit = lambda track: track.stdT/track.stdL > .8 and track.nSavedPix > 9 and track.errT/track.L < 1.5


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
                        xKey = np.array( [ x[j] for j in range(len(self.tracks)) if guess2(j,key) ] )
                        xTrain = np.array( [ x[j] for j in values ] )
                        m = np.mean(xKey)
                        ax_guess.hist( xKey, bins = bins, log=True, histtype = 'step', color = colormapping[key], label = key )
                        ax_train.hist( xTrain, bins = bins, log=True, histtype = 'step', color = colormapping[key], label = key )
                    ax_guess.legend()
                    
                    #ax.set_ylabel(var2)
                else:
                    print '\t', i2, var2
                    y = np.array(self.__dict__[var2])
                    ax_guess.scatter( y, x, s = s, c = [ colormapping[key] for i in range(len(self.tracks)) for key in categories.keys() if key != 'u' and guess2(i,key) ], linewidths = 0 )
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
####################################################################################################################
class Track:
####################################################################################################################
    def __init__( self, **kwargs ):
        for name, item in kwargs.items():
            self.__dict__[name] = item
####################################################################################################################
    def check( self ):
        self.xPix = np.array(self.xPix)
        self.yPix = np.array(self.yPix)
        self.ePix = np.array(self.ePix)
        self.ePix[ self.ePix < 0 ] = 0
        
        return len(self.xPix)
####################################################################################################################
    def eigcov( self ):
        self.check()
        #C =  np.cov( self.xPix, self.yPix )
        C = cov( np.vstack((self.xPix, self.yPix)), self.ePix )
        
        w, v = np.linalg.eig( C )
        
        if w[0] > w[1]: self.wLong, self.wTrans, self.vLong, self.vTrans = w[1], w[0], v[:,1], v[:,0]
        else: self.wLong, self.wTrans, self.vLong, self.vTrans = w[0], w[1], v[:,0], v[:,1]
        
        self.lPix = np.inner( self.vLong, zip( self.xPix, self.yPix ) )
        self.L = max(self.lPix) - min(self.lPix)
        self.lReduced = (self.lPix - min(self.lPix))/(max(self.lPix) - min(self.lPix))
        self.tPix = np.inner( self.vTrans, zip( self.xPix, self.yPix ) )
        self.T = max(self.tPix) - min(self.tPix)
        self.stdL = std( self.lPix, weights = self.ePix )
        self.stdT = std( self.tPix, weights = self.ePix )

        tMean = np.average( self.tPix, weights = self.ePix )
        lMean = np.average( self.lPix, weights = self.ePix )
        if np.max(self.ePix) == 0:
            self.errT = 0
            self.errL = 0
        else:
            self.errT = np.sum( ( np.exp( -.5* (self.tPix - tMean)**2/self.stdT**2) - self.ePix/np.max(self.ePix) )**2 )
            self.errL = np.sum( ( np.exp( -.5* (self.lPix - lMean)**2/self.stdL**2 ) - self.ePix/np.max(self.ePix) )**2 )

        #if self.stdL < self.stdT: 
            #self.stdL, self.stdT = self.stdT, self.stdL
            #self.L, self.T = self.T, self.L
            #self.errL, self.errT = self.errT, self.errL
        
        self.xCenter = np.average( self.xPix )
        self.yCenter = np.average( self.yPix )
        
        self.r2 = (self.xPix - self.xCenter)**2 + (self.yPix - self.yCenter)**2
        self.phi = np.arctan2( self.yPix - self.yCenter, self.xPix - self.xCenter )
        self.normalization = np.sum( self.r2*self.ePix )
        self.eccentricity = {}
        self.eccentricityRaw = {}
        self.eccentricityRawU = {}
        self.computeEccentricity(2)
        self.computeEccentricity(3)
        self.computeEccentricity(4)
        self.ecc2 = self.eccentricity[2]
        self.ecc3 = self.eccentricity[3]
        self.ecc4 = self.eccentricity[4]
        self.computeSigmasOpt2(10)
        self.stdZ10, self.stdZ010, self.stdZerr10 = self.stdZ, self.stdZ0, self.stdZerr
        self.stdZ_ex10, self.stdZ0_ex10, self.stdZerr_ex10 = self.stdZ_ex, self.stdZ0_ex, self.stdZerr_ex
        self.slope10, self.intp10, self.slopeErr10 = self.slope, self.intp, self.slopeErr
        self.computeSigmasOpt2(30)
        self.stdZ30, self.stdZ030, self.stdZerr30 = self.stdZ, self.stdZ0, self.stdZerr
        self.stdZ_ex30, self.stdZ0_ex30, self.stdZerr_ex30 = self.stdZ_ex, self.stdZ0_ex, self.stdZerr_ex
        self.slope30, self.intp30, self.slopeErr30 = self.slope, self.intp, self.slopeErr
        self.computeSigmasOpt2(50)
        self.stdZ50, self.stdZ050, self.stdZerr50 = self.stdZ, self.stdZ0, self.stdZerr
        self.stdZ_ex50, self.stdZ0_ex50, self.stdZerr_ex50 = self.stdZ_ex, self.stdZ0_ex, self.stdZerr_ex
        self.slope50, self.intp50, self.slopeErr50 = self.slope, self.intp, self.slopeErr
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
        mat = ax.matshow( self.draw( np.zeros( (max(self.xPix)-min(self.xPix)+1, max(self.yPix)-min(self.yPix)+1) ), log = log ), cmap = plt.get_cmap(cmap), **kwargs )
        #plt.colorbar(mat, fraction = 0.046, pad = 0.04)
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
    print 'catalog file', fname
    opt = False
    catalog = Catalog(fname, opt = opt )
    if not opt: catalog.write()
    #matrix = np.zeros( (2048, 4096) )

    #t1 = timer('eigcov')
    #for track in catalog.getTracks():
        #track.eigcov()
        #track.computeEccentricity(2)
        #track.computeEccentricity(3)
        #track.computeEccentricity(4)
    #t1.factor( text = 'per track', f = 1./len(catalog.getTracks()) )
    #del(t1)
    #catalog.write()

if __name__ == "__main__":
    if '--update' in sys.argv:
        print "update"
        mainUpdate( sys.argv[-1] )
    elif '--split' in sys.argv:
        print "split"
        catalog = Catalog2()
        catalog.listRunID( sys.argv[-1] )
    elif '--analyze' in sys.argv:
        print "analyze"
        catalog = Catalog(sys.argv[-1], complete = True)
        catalog.PCA()
    elif '--learn' in sys.argv:
        print "learn"
        catalog = Catalog(sys.argv[-1], complete = True)
        catalog.learn()
    else:
        print 'no recognized option given'
        mainHelp()
    exit(0)


def plot_transverseLongHistogram( catalog ):
    f, ax = plt.subplots(1)
    ax.set_title( 'tranverse vs longitudinal histogram' )
    ax.set_xlabel( 'trans long' )
    ax.set_ylabel( 'count' )
    num_bins = 100
    ax.hist( [ track.eigRatio() for track in catalog.getTracks() ], num_bins, log = 'y', histtype='step' )
    ax.hist( [ track.sigmaRatio() for track in catalog.getTracks() ], num_bins, log = 'y', histtype='step' )
    plt.savefig('ratio_hist.png')
    plt.clf()
    f, ax = plt.subplots(1)
    ax.set_title( 'longitudinal histogram' )
    ax.set_xlabel( 'long' )
    ax.set_ylabel( 'count' )
    num_bins = 100
    ax.hist( [ np.sqrt(track.wLong) for track in catalog.getTracks() ], num_bins, log = 'y', histtype='step' )
    ax.hist( [ track.stdL for track in catalog.getTracks() ], num_bins, log = 'y', histtype='step' )
    plt.savefig('long_hist.png')
    plt.clf()
    f, ax = plt.subplots(1)
    ax.set_title( 'tranverse histogram' )
    ax.set_xlabel( 'trans' )
    ax.set_ylabel( 'count' )
    num_bins = 100
    n, bins, pathces = ax.hist( [ np.sqrt( track.wTrans ) for track in catalog.getTracks() ], num_bins, log = 'y', histtype='step' )
    n, bins, pathces = ax.hist( [ track.stdT for track in catalog.getTracks() ], num_bins, log = 'y', histtype='step' )
    #order = np.argsort( n )
    #print n[order[-1]], bins[order[-1]]
    plt.savefig('trans_hist.png')
    plt.clf()

def plot_transverseLongScatter( catalog ):
    f, ax = plt.subplots(1)
    ax.set_title( 'tranverse vs longitudinal eigenvalues' )
    ax.set_xlabel( 'trans/long' )
    ax.set_ylabel( 'log(trans)' )
    x = np.array([ np.sqrt(track.wLong) for track in catalog.getTracks() ])
    y = np.array([ np.sqrt(track.wTrans) for track in catalog.getTracks() ])
    plot = ax.plot
    plot( np.log(y/x), np.log(y), 'k.' )
    plt.savefig('scatter_logtransverse_vs_ratio.png')
    plt.clf()
    f, ax = plt.subplots(1)
    ax.set_title( 'tranverse vs longitudinal eigenvalues' )
    ax.set_xlabel( 'long' )
    ax.set_ylabel( 'log(trans)' )
    #x2 = np.r_[ 0 : np.max(y) : 100j ]
    plot = ax.plot
    plot( np.array(x), np.log(y), 'k.' )
    plt.savefig('scatter_logtransverse_vs_long.png')
    plt.clf()

    f, ax = plt.subplots(1)
    ax.set_title( 'tranverse vs longitudinal std' )
    ax.set_xlabel( 'trans/long' )
    ax.set_ylabel( 'log(trans)' )
    x = np.array([ track.stdL for track in catalog.getTracks() ])
    y = np.array([ track.stdT for track in catalog.getTracks() ])
    plot = ax.plot
    plot( np.log(y/x), np.log(y), 'k.' )
    plt.savefig('scatter_logtransverse_vs_ratio_std.png')
    plt.clf()
    f, ax = plt.subplots(1)
    ax.set_title( 'tranverse vs longitudinal std' )
    ax.set_xlabel( 'long' )
    ax.set_ylabel( 'trans' )
    plot = ax.plot
    plot( np.array(x), np.log(y), 'k.' )
    plt.savefig('scatter_transverse_vs_long_std.png')
    plt.clf()

def plot_transverseLongHist2d( catalog ):
    f, ax = plt.subplots(1)
    ax.set_title( 'tranverse vs longitudinal eigenvalues' )
    ax.set_xlabel( 'trans/long' )
    ax.set_ylabel( 'log(trans)' )
    #ax.set_clabel( 'count' )
    x = np.array([ np.sqrt(track.wLong) for track in catalog.getTracks() ])
    y = np.array([ np.sqrt(track.wTrans) for track in catalog.getTracks() ])
    ax.hist2d( y/x, np.log(y), bins = len(catalog.getTracks())**(1./3), cmap = 'rainbow', norm=matplotlib.colors.LogNorm() )
    #plt.colorbar()
    plt.savefig('hist2d_transverse_vs_ratio.png')
    plt.clf()
    
    f, ax = plt.subplots(1)
    ax.set_title( 'tranverse vs longitudinal std' )
    ax.set_xlabel( 'trans/long' )
    ax.set_ylabel( 'log(trans)' )
    #ax.set_clabel( 'count' )
    x = np.array([ track.stdL for track in catalog.getTracks() ])
    y = np.array([ track.stdT for track in catalog.getTracks() ])
    ax.hist2d( y/x, np.log(y), bins = len(catalog.getTracks())**(1./3), cmap = 'rainbow', norm=matplotlib.colors.LogNorm() )
    #plt.colorbar()
    plt.savefig('hist2d_transverse_vs_long_std.png')
    plt.clf()
    
#plot_transverseLongScatter( catalog )
plot_transverseLongHist2d( catalog )
#plot_transverseLongHistogram( catalog )

exit(0)

def computeAlpha( tracks, nbins ):
    results = [ track.computeSigmas( nbins ) for track in tracks ]
    stdTBin, bins, slopes, intercepts, std_errs = map( np.array, zip(*results) )
    weights = np.zeros_like(std_errs)
    weights[ std_errs > 0 ] = 1./std_errs[ std_errs > 0]
    alphaN, betaN, tmp, tmp2, std_errN = stats.linregress( bins[0,:], np.average( stdTBin, axis = 0 ) )

    resultsOpt = [ track.computeSigmasOpt( min_pixels = 50 ) for track in tracks ]
    stdTBinOpt, binsOpt, slopesOpt, interceptsOpt, std_errsOpt = map( np.array, zip(*resultsOpt) )

    resultsOpt2 = [ track.computeSigmasOpt( min_pixels = 20 ) for track in tracks ]
    stdTBinOpt2, binsOpt2, slopesOpt2, interceptsOpt2, std_errsOpt2 = map( np.array, zip(*resultsOpt2) )
    
    #for track in tracks: track.computeSigmaTrans( z = zAxis ).computeSlope( z = zAxis )
    #slopes = np.array( [ np.abs(track.slope0) for track in tracks ] )
    #slopes1 = np.array( [ np.abs(track.slope1) for track in tracks ] )
    #std_errs = np.array( [ track.std_err0 for track in tracks ] )
    #std_errs1 = np.array( [ track.std_err1 for track in tracks ] )
    
    #stdTZ = np.array( [ track.stdTZ for track in tracks ])
    #stdTZ1 = np.array( [ track.stdTZ1 for track in tracks ])
    
    #sigmaAve = np.average( stdTZ, axis = 0 )
    #sigmaErr = stats.sem( stdTZ, axis = 0 )
    ##w = 1./sigmaErr
    #alphaN, betaN, tmp, tmp2, std_errN = stats.linregress( zAxis, sigmaAve )
    ##alpha, beta, tmp, tmp2, std_err = stats.linregress( w*zAxis, w*sigmaAve )

    #sigmaAve1 = np.average( stdTZ1, axis = 0 )
    #sigmaErr1 = stats.sem( stdTZ1, axis = 0 )
    ##w = 1./sigmaErr1
    #alpha1N, beta1N, tmp, tmp2, std_err1N = stats.linregress( zAxis, sigmaAve1 )
    ##alpha1, beta1, tmp, tmp2, std_err1 = stats.linregress( w*zAxis, w*sigmaAve1 )

    #w = 1./std_errs
    #w1 = 1./std_errs1
    return {'slope_ave_N': [alphaN, betaN, std_errN ], 
            'ave_slope': [ np.average(slopes, weights = weights ), np.average(intercepts, weights = weights ), std( slopes, weights = weights )/np.sqrt(len(slopes)) ],
            'ave_slope_N': [ np.average( slopes ), np.average( intercepts ), np.std( slopes )/np.sqrt(len(slopes)) ],
            'slopes': [slopes, intercepts, std_errs],
            'slopesOpt': [slopesOpt, interceptsOpt, std_errsOpt],
            'ave_slopeOpt_N': [ np.average( slopesOpt ), np.average( interceptsOpt ), np.std( slopesOpt )/np.sqrt(len(slopesOpt)) ],
            'slopesOpt2': [slopesOpt2, interceptsOpt2, std_errsOpt2],
            'ave_slopeOpt2_N': [ np.average( slopesOpt2 ), np.average( interceptsOpt2 ), np.std( slopesOpt2 )/np.sqrt(len(slopesOpt2)) ],
            #'z':  bin,
            'N': len(slopes),
            'tracks': [bins, stdTBin ], 
            #'excl_slope_ave_N': [alpha1N, std_err1N],
            #'excl_ave_slope': [np.average( w1**2*slopes1 )/np.average(w1**2), stats.sem( w1*slopes1 )/np.average(w1**2)],
            #'excl_ave_slope_N': [np.average( slopes1 ), stats.sem( slopes1 )],
            #'excl_ave_tracks': [sigmaAve1, sigmaErr1]
            }

def scan_stdRatios( catalog, zAxis ):
    step = .02
    stdRatios = np.arange( 0.07, 1, step )
    
    result = np.array( [  computeAlpha( [ track for track in catalog.getTracks() if track.sigmaRatio() < stdRatio + step/2 and track.sigmaRatio() > stdRatio - step/2 ], zAxis  ) 
                        for stdRatio in stdRatios ] )
    resultCum = np.array( [  computeAlpha( [ track for track in catalog.getTracks() if track.sigmaRatio() < stdRatio ], zAxis  ) 
                        for stdRatio in stdRatios ] )
    
    f, ax = plt.subplots(1)
    ax.set_title( 'witdh/depth correlation per tranverse/longitudinal' )
    ax.set_xlabel( 'trans/long' )
    ax.set_ylabel( 'width/depth' )
    ax.errorbar( stdRatios, result[:,0], xerr = step/2, yerr = result[:,1], fmt='.' )
    ax.errorbar( stdRatios, result[:,2], xerr = step/2, yerr = result[:,3], fmt='.' )
    plt.savefig('width_vs_depth.png')
    plt.clf()
    f, ax = plt.subplots(1)
    ax.set_title( 'witdh/depth correlation per tranverse/longitudinal cumulative' )
    ax.set_xlabel( 'trans/long' )
    ax.set_ylabel( 'width/depth' )
    ax.errorbar( stdRatios, resultCum[:,0], yerr = resultCum[:,1], fmt='.' )
    ax.errorbar( stdRatios, resultCum[:,2], yerr = resultCum[:,3], fmt='.' )
    plt.savefig('width_vs_depthCum.png')
    plt.clf()

def scan_Bins( catalog, bins_, stdRatio, TransCut = 100, plot = False, plot2 = False ):
    
    if False:
        resultAll = [ computeAlpha( catalog.getTracks(), nbin ) for nbin in bins_ ]
        for nbin, res in zip( bins_, resultAll ):
            f, ax = plt.subplots(1)
            ax.set_title( 'slope histogram %s' % nbin )
            ax.set_ylabel( 'count' )
            ax.set_xlabel( 'slope' )
            #
            slopes = np.array(res['slopes'][0])
            slopes[slopes!=slopes] = 0
            ax.hist( slopes, np.sqrt(len(catalog.getTracks())), log = 'y', histtype = 'step' )
            plt.savefig('hist_slope%s.png' % nbin )
            plt.clf()
            plt.close(f)
            f, ax = plt.subplots(1)
            ax.set_title( 'slope versus stdRatio %s' % nbin )
            ax.set_ylabel( 'slope' )
            ax.set_xlabel( 'stdRatio' )
            ax.scatter( [ track.sigmaRatio() for track in catalog.getTracks() ], res['slopes'][0] )
            plt.savefig('scatter_slope_k%s.png' % nbin )
            plt.clf()
            plt.close(f)
            f, ax = plt.subplots(1)
            ax.set_title( 'slope versus stdTrans %s' % nbin )
            ax.set_ylabel( 'slope' )
            ax.set_xlabel( 'stdTrans' )
            ax.scatter( [ track.stdT for track in catalog.getTracks() ], res['slopes'][0] )
            plt.savefig('scatter_slope_trans%s.png' % nbin )
            plt.clf()
            plt.close(f)
    
    selectedTracks = [ track for track in catalog.getTracks() if track.sigmaRatio() < stdRatio and track.stdT < TransCut ]
    result = [ computeAlpha( selectedTracks, nbin ) for nbin in bins_ ]
    if plot:
        for i,track in enumerate(selectedTracks):
            f, ax = plt.subplots(1)
            ax.set_title( r'$k=%.3f$, $\sigma_\perp=%.3f$; ratio$=%.3f$, long$=%.3f$'%( track.sigmaRatio(), track.stdT, track.eigRatio(), np.sqrt(track.wTrans) ) )
            track.plot(ax=ax, log = True, cmap = 'rainbow')
            Lx = np.amax(track.xPix) - np.amin(track.xPix)
            Ly = np.amax(track.yPix) - np.amin(track.yPix)
            #ax.arrow( Ly/2, Lx/2, 20*track.wTrans[1], 20*track.wTrans[0] )
            #ax.arrow( Ly/2, Lx/2, 20*track.wLong[1], 20*track.wLong[0] )
            
            nbin = 10
            for n in np.linspace( 0.1, 0.9, nbin+1 ):
                nn = 2*n-1
                dx = Lx/2
                dy = Ly/2
                ax.arrow( Ly/2 + nn*dy*track.vLong[1], Lx/2 + nn*dx*track.vLong[0], 100*track.vTrans[1], 100*track.vTrans[0] )
                ax.arrow( Ly/2 + nn*dy*track.vLong[1], Lx/2 + nn*dx*track.vLong[0], -100*track.vTrans[1], -100*track.vTrans[0] )
            
                
            #ax[1].hist2d(track.xPix, track.yPix, bins = [max(track.xPix)-min(track.xPix)+1,max(track.yPix)-min(track.yPix)+1], weights = track.ePix)
            #plt.colorbar()
            plt.savefig('track%03i.png'%i)
            plt.clf()
            plt.close(f)
            
            f, ax = plt.subplots(len(bins_))
            ax[0].set_title( r'$k=%.3f$, $\sigma_\perp=%.3f$; ratio$=%.3f$, long$=%.3f$'%( track.sigmaRatio(), track.stdT, track.eigRatio(), np.sqrt(track.wTrans) ) )
            for n,nbin in enumerate(bins_):
                ax[n].plot( result[n]['tracks'][0][i], result[n]['tracks'][1][i], 'k.' )
                ax[n].plot( result[n]['tracks'][0][i], result[n]['slopes'][0][i]*result[n]['tracks'][0][i]+ result[n]['slopes'][1][i], 'r-', label = r'$%.3f,%.3f$' %(result[n]['slopes'][0][i],result[n]['slopes'][1][i]) )
                ax[n].plot( result[n]['tracks'][0][i], result[n]['slopesOpt'][0][i]*result[n]['tracks'][0][i]+ result[n]['slopesOpt'][1][i], 'r-', label = r'$%.3f,%.3f$' %(result[n]['slopesOpt'][0][i],result[n]['slopesOpt'][1][i]) )
                ax[n].legend()
            plt.savefig('track%03i_slope.png'%i)
            plt.clf()
            plt.close(f)
            


    f, ax = plt.subplots(1)
    ax.set_title( 'slope per number of bins k=%s t=%s N%s' % (stdRatio, TransCut, len(selectedTracks)) )
    ax.set_xlabel( 'bins' )
    ax.set_ylabel( 'slope' )
    #ax.errorbar( bins_, [ res['slope_ave_N'][0] for res in result ], yerr = [ res['slope_ave_N'][2] for res in result ], fmt='.', label='slope ave N' )
    ax.errorbar( bins_, [ res['ave_slope'][0] for res in result ], yerr = [ res['ave_slope'][2] for res in result ], fmt='.', label='ave slope' )
    ax.errorbar( bins_, [ res['ave_slope_N'][0] for res in result ], yerr = [ res['ave_slope_N'][2] for res in result ], fmt='.', label='ave slope N' )
    ax.errorbar( bins_, [ res['ave_slopeOpt_N'][0] for res in result ], yerr = [ res['ave_slopeOpt_N'][2] for res in result ], fmt='.', label='ave slope opt N' )
    ax.errorbar( bins_, [ res['ave_slopeOpt2_N'][0] for res in result ], yerr = [ res['ave_slopeOpt2_N'][2] for res in result ], fmt='.', label='ave slope opt2 N' )
    #ax.errorbar( bins_, [ res['excl_slope_ave'][0] for res in result ], yerr = [ res['excl_slope_ave'][1] for res in result ], fmt='.', label='excl slope ave' )
    #ax.errorbar( bins_, [ res['excl_slope_ave_N'][0] for res in result ], yerr = [ res['excl_slope_ave_N'][1] for res in result ], fmt='.', label='excl slope ave N' )
    #ax.errorbar( bins_, [ res['excl_ave_slope'][0] for res in result ], yerr = [ res['excl_ave_slope'][1] for res in result ], fmt='.', label='excl ave slope' )
    #ax.errorbar( bins_, [ res['excl_ave_slope_N'][0] for res in result ], yerr = [ res['excl_ave_slope_N'][1] for res in result ], fmt='.', label='excl ave slope N' )
    ax.legend()# loc=2 )
    plt.savefig('width_vs_depth_bink%sP%sN%s.png' % (stdRatio,TransCut, len(selectedTracks)) )
    plt.close(f)
    
    f, ax = plt.subplots(1)
    ax.set_title( 'hist slope k=%s t=%s N%s' % (stdRatio, TransCut, len(selectedTracks)) )
    ax.set_xlabel( 'slope' )
    ax.set_ylabel( 'count' )
    ax.hist( result[0]['slopesOpt'][0], bins = np.sqrt(len(selectedTracks)), normed = True, histtype = 'step', label='opt' )
    ax.hist( result[0]['slopesOpt2'][0], bins = np.sqrt(len(selectedTracks)), normed = True, histtype = 'step', label='opt' )
    for nbins, res in zip(bins_, result):
        ax.hist( res['slopes'][0], bins = np.sqrt(len(selectedTracks)), normed = True, histtype = 'step', label='%s'%nbins )
    #ax.errorbar( bins_, [ res['excl_slope_ave'][0] for res in result ], yerr = [ res['excl_slope_ave'][1] for res in result ], fmt='.', label='excl slope ave' )
    #ax.errorbar( bins_, [ res['excl_slope_ave_N'][0] for res in result ], yerr = [ res['excl_slope_ave_N'][1] for res in result ], fmt='.', label='excl slope ave N' )
    #ax.errorbar( bins_, [ res['excl_ave_slope'][0] for res in result ], yerr = [ res['excl_ave_slope'][1] for res in result ], fmt='.', label='excl ave slope' )
    #ax.errorbar( bins_, [ res['excl_ave_slope_N'][0] for res in result ], yerr = [ res['excl_ave_slope_N'][1] for res in result ], fmt='.', label='excl ave slope N' )
    ax.legend()# loc=2 )
    plt.savefig('hist_slopes_bink%sP%sN%s.png' % (stdRatio,TransCut, len(selectedTracks)) )
    plt.close(f)
    
    alpha, beta = np.average( [result[0]['slopesOpt'][0], result[0]['slopesOpt2'][0]] ), np.average( [result[0]['slopesOpt'][1], result[0]['slopesOpt2'][1]] )
    print 'alpha', alpha, beta
    selectedTracks2 = [ track for track in catalog.getTracks() if track.sigmaRatio() > .9 ] #and track.stdT < 1 ]
    
    if plot2:
        for i,track in enumerate(selectedTracks2):
            f, ax = plt.subplots(1)
            ax.set_title( r'$k=%.3f$, $\sigma_\perp=%.3f$; ratio$=%.3f$, long$=%.3f$'%( track.sigmaRatio(), track.stdT, track.eigRatio(), np.sqrt(track.wTrans) ) )
            track.plot(ax=ax, log = True, cmap = 'rainbow')
            plt.savefig('track2_%03i.png'%i)
            plt.clf()
            plt.close(f)

    z = np.array( [ (track.stdT - beta)/alpha for track in selectedTracks2 ] )
    f, ax = plt.subplots(1)
    ax.set_title( 'z k=%s t=%s N%s' % (stdRatio, TransCut, len(selectedTracks)) )
    ax.set_xlabel( 'z' )
    ax.set_ylabel( 'count' )
    ax.hist( z, bins = np.sqrt(len(selectedTracks2)), normed = True, histtype = 'step' )
    plt.savefig('hist_zk%sP%sN%s.png' % (stdRatio,TransCut, len(selectedTracks)) )
    plt.close(f)
    

#exit(0)
nbins = [ 5, 10, 25, 50, 75]#, 100, 150 ]
#nbins = [ 25 ] #, 50, 75, 100, 150 ]
#for k in [0.07,0.1,.2,.3]:
    #scan_Bins( catalog, nbins, k )
#scan_Bins( catalog, nbins, stdRatio = .04, plot = True )
scan_Bins( catalog, nbins, stdRatio = .05, plot = False, TransCut = 1 )
scan_Bins( catalog, nbins, stdRatio = .04, plot = False, TransCut = .7 )
#scan_Bins( catalog, nbins, stdRatio = .1 )
#scan_Bins( catalog, nbins, stdRatio = .2 )
#scan_Bins( catalog, nbins, .1, TransCut = 3 )
#scan_Bins( catalog, nbins, .2, TransCut = 3 )
#scan_Bins( catalog, nbins, .3, TransCut = 3 )
#scan_Bins( catalog, nbins, .1, TransCut = 2.8 )
#scan_Bins( catalog, nbins, .2, TransCut = 2.8 )
#scan_Bins( catalog, nbins, .3, TransCut = 2.8 )
#scan_Bins( catalog, nbins, .3, TransCut = 2.6 )
#scan_Bins( catalog, nbins, .2, TransCut = 2.6 )
#scan_Bins( catalog, nbins, .1, TransCut = 2.6 )
#zAxis = np.linspace( 0.1, 0.9, 10+1 )
#scan_stdRatios( catalog, zAxis )

#def compute

exit(0)



