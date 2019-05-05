import sys
import root_numpy
import numpy as np
import scipy.stats
import scipy.optimize
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Circle

usagestr="""\
    usage:
    h       this message
    f       input root file
    x       expression for xAxis
    logy    sets log scale
"""

def cov(m, a):
    v1 = np.sum(a) # N
    m -= np.sum(m * a, axis=1, keepdims=True) / v1 # m = m - (sum m)/N
    return np.dot(m * a, m.T) / v1 # m*m.T * N/(N**2 - (N-1)*N)

def circularMask(x, y, centerx, centery, radius):
    return (x-centerx)**2 + (y-centery)**2 < radius**2

def std( x, weights=None ):
    """
    calculate the wieghted standard deviation by using the numpy average function
    """
    x=np.array(x)
    if weights is None: weights = np.ones_like(x)
    weights[weights!=weights] = 0
    std2 = np.average( x**2, weights = weights ) - np.average( x, weights = weights )**2
    return np.sqrt( np.abs(std2) )

def fitMuon( fname, index ):
    data = root_numpy.root2array(fname, treename = 'hitSumm', start=index, stop=index+1 )
    xPix = data['xPix'][0].astype(float)
    yPix = data['yPix'][0].astype(float)
    ePix = data['ePix'][0].astype(float)
    ePix -= np.min(ePix)
    
    print 'len', len(ePix)
    C = cov( np.vstack( (xPix, yPix) ), ePix )
    w_, v_ = np.linalg.eigh(C)
    print v_[:,1]
    ediffPix = np.zeros_like(ePix)

    tPix = np.dot( zip(xPix,yPix), v_[:,0] )
    lPix = np.dot( zip(xPix,yPix), v_[:,1] )
    tPix -= np.min(tPix)
    lPix -= np.min(lPix)
    
    meanT = np.average( tPix, weights=ePix )
    meanL = np.average( lPix, weights=ePix )
    stdT = std( tPix, weights=ePix )
    L = np.max(lPix)
    print 'props', meanT, meanL, stdT, L
    p12 = (0, 1, 0, 1)
    func2 = lambda x,a,b,c,d: scipy.stats.norm.pdf( x[0]-meanT, a, b*stdT + c*(x[1]-meanL)/L)*d*np.sum(ePix)
    try: pp2 = scipy.optimize.curve_fit( func2, np.vstack((tPix,lPix)), ePix, p0=p12 )
    except: return 0, 0, L
    if pp2[0][2] != pp2[0][2]: return 0, 0, L
    print 'first estimate', pp2[0]
    ediffPix = (ePix-func2(np.vstack((tPix,lPix)), *pp2[0]))**2
    mean_ediff, std_ediff = np.mean(ediffPix), np.std(ediffPix)
    print mean_ediff, std_ediff
    estCount = 1
    mask = np.abs(ediffPix-mean_ediff)<std_ediff
    looptPix, looplPix, loopePix = tPix[mask], lPix[mask], ePix[mask]
    err = 1
    while 1:
        estCount += 1
        pp22 = pp2[0][2]
        try: pp2 = scipy.optimize.curve_fit( func2, np.vstack((looptPix,looplPix)), loopePix, p0=pp2[0] )
        except: return err, 0, L
        if pp2[0][2] != pp2[0][2]: return err, pp22, L
        print estCount,'estimate', pp2[0]
        #if pp2[0][2] == 0: break
        ediffPix = (loopePix-func2(np.vstack((looptPix,looplPix)), *pp2[0]))**2
        mean_ediff0, std_ediff0 = mean_ediff, std_ediff
        mean_ediff, std_ediff = np.mean(ediffPix), np.std(ediffPix)
        #print 'pp2', pp2[0]
        err = np.abs(pp22-pp2[0][2])/np.abs(pp22)
        print 'err', err, len(loopePix), len(ePix)
        if err < 0.01: 
            break
        if estCount > 10:
            print 'did not converge'
            break
        mask = np.abs(ediffPix-mean_ediff)<std_ediff
        looptPix, looplPix, loopePix = looptPix[mask], looplPix[mask], loopePix[mask]
    return err, pp2[0][2], L

def scanTracks( fname, minindex, maxindex):
    ret = [ fitMuon(fname, index) for index in range(minindex,maxindex) ]
    print ret
    fig = plt.figure(figsize=(10,10))
    fig.subplots_adjust(hspace=0, wspace=0)
    ax = fig.add_subplot(111)
    print 
    ret2 = zip(*ret)
    y = np.abs(ret2[1])
    y[y>2]=0
    ax.plot( ret2[2], y, '.' )
    fig.subplots_adjust()
    fig.savefig('viewer.pdf')
    return

def muonFit( x, *p ):
    y = np.zeros(len(x.T))
    meanx = (p[0]+p[1])/2
    diffx = (p[0]-p[1])/2
    mask1 = np.abs(x[1]-meanx)<np.abs(diffx)
    y[ mask1 ] = p[2]*scipy.stats.norm.pdf(x[0][mask1]-p[3], 0, p[4] + p[5]*x[1][mask1] )
    #y = p[2]*scipy.stats.norm.pdf(x[0]-p[3], 0, p[4] + p[5]*x[1] )
    mask2 = x[1] < p[0]
    y[ mask2 ] = p[2]*scipy.stats.norm.pdf(np.sqrt( (x[0][mask2]-p[3])**2 + (x[1][mask2]-p[0])**2 ), 0, p[4] + p[5]*p[0] )
    mask3 = x[1] > p[1]
    y[ mask3 ] = p[2]*scipy.stats.norm.pdf(np.sqrt( (x[0][mask3]-p[3])**2 + (x[1][mask3]-p[1])**2 ), 0, p[4] + p[5]*p[1] )
    
    return y

def plotTrack( fname, index ):
    print 'plotTrack', index
    baseSel = 'flag==0' + ('&&%s<%s'%(expr,Max) if Max else '') + ('&&%s>%s'%(expr,Min) if Min else '')
    data = root_numpy.root2array(fname, treename = 'hitSumm', start=index, stop=index+1 )
    
    divisions=4
    
    fig = plt.figure(figsize=(10,10))
    fig.subplots_adjust(hspace=0, wspace=0)
    fig2 = plt.figure(figsize=(10,10))
    fig2.subplots_adjust(hspace=0, wspace=0)
    fig3 = plt.figure(figsize=(10,10))
    fig3.subplots_adjust(hspace=0, wspace=0)

    #xPix = np.array(data['xPix'][0], dtype='float32')
    #yPix = np.array(data['yPix'][0], dtype='float32')
    #ePix = np.array(data['ePix'][0], dtype='float32')
    xPix = data['xPix'][0].astype(float)
    yPix = data['yPix'][0].astype(float)
    ePix = data['ePix'][0].astype(float)
    ePix -= np.min(ePix)
    
    offset = (np.min(data['xPix'][0]), np.min(data['yPix'][0]))
    matrix = np.zeros( (np.max(data['xPix'][0])-np.min(data['xPix'][0])+1, np.max(data['yPix'][0])-np.min(data['yPix'][0])+1 ) )
    #for x,y,e in zip( data['xPix'][0], data['yPix'][0], data['ePix'][0] ): matrix[x-offset[0], y-offset[1]] += e#np.log10(e+1) if e>0 else 0
    matrix[ (xPix - np.min(xPix)).astype(int), (yPix - np.min(yPix)).astype(int) ] = ePix
    gradmatrixx = (matrix[:-2,1:-1]-matrix[2:,1:-1])
    gradmatrixy = (matrix[1:-1,:-2]-matrix[1:-1,2:])
        
    gradd = (matrix[0:-2,0:-2]-matrix[2:,2:])/np.sqrt(2)
    gradD = (matrix[0:-2,2:]-matrix[2:,0:-2])/np.sqrt(2)
    
    gradmatrixx = (gradmatrixx + gradd + gradD)/2
    gradmatrixy = (gradmatrixy + gradd - gradD)/2
    
    ax = fig.add_subplot(311)
    ax.imshow( matrix.T, cmap=plt.cm.YlOrRd, origin='lower' )

    #C = cov( np.vstack( (data['xPix'][0].astype(float), data['yPix'][0].astype(float)) ), data['ePix'][0] )
    C = cov( np.vstack( (xPix, yPix) ), ePix )
    w_, v_ = np.linalg.eigh(C)
    #print C.shape
    #print w_.shape, w_
    print v_[:,1]
    angle0 = np.angle(v_[0,1]+v_[1,1]*1.j, deg=True)
    angle0rad = np.angle(v_[0,1]+v_[1,1]*1.j)
    #print 'angle0', angle0

    flow = lambda n,phi,x,y: np.sum( np.exp(1.j*n*(np.angle(x+y*1.j)-phi) )*np.absolute(x+y*1.j) )
    v0 = flow(0,0,gradmatrixx,gradmatrixy)
    #print flow(0,0,gradmatrixx,gradmatrixy), flow(1,0,gradmatrixx,gradmatrixy), flow(2,0,gradmatrixx,gradmatrixy)
    psi1 = np.angle(flow(1,0,gradmatrixx,gradmatrixy))
    psi2 = np.angle(flow(2,0,gradmatrixx,gradmatrixy))/2
    #print np.angle(flow(1,0,gradmatrixx,gradmatrixy)), np.angle(flow(2,0,gradmatrixx,gradmatrixy)), angle0rad
    #print flow(2,psi2,gradmatrixx,gradmatrixy)
    #print np.absolute(flow(1,0,gradmatrixx,gradmatrixy)/v0 ), np.absolute( flow(2,0,gradmatrixx,gradmatrixy)/v0), np.absolute( flow(3,0,gradmatrixx,gradmatrixy)/v0 ), np.absolute( flow(4,0,gradmatrixx,gradmatrixy)/v0 ), np.absolute( flow(5,0,gradmatrixx,gradmatrixy)/v0 ), np.absolute( flow(6,0,gradmatrixx,gradmatrixy)/v0 )

    cx, cy = matrix.shape[0]/2, matrix.shape[1]/2
    #ax.quiver( [cx], [cy], [w_[1]*v_[0,1]], [w_[1]*v_[1,1] ], units='dots', scale=.5, color='r' )
    #ax.quiver( [matrix.shape[0]/2], [matrix.shape[1]/2], [w_[1]*v_[0,1]], [w_[1]*v_[1,1] ], units='dots', scale=.5, color='r' )
    
    #ax.quiver( [cx], [cy], [np.cos(psi1)], [np.sin(psi1)], units='dots', scale=.05, color='b' )
    #ax.quiver( [cx,cx], [cy,cy], [np.cos(psi2), np.cos(psi2+np.pi)], [np.sin(psi2), np.sin(psi2+np.pi)], units='dots', scale=.05, color='g' )

    X, Y = np.mgrid[1:gradmatrixx.shape[0]+1, 1:gradmatrixx.shape[1]+1]
    ax.quiver( X, Y, gradmatrixx, gradmatrixy, units='dots', scale=500 )

    axh = fig.add_subplot(312)
    angles = np.angle( gradmatrixx.flatten()+gradmatrixy.flatten()*1.j, deg=True )
    angles = angles - psi2/np.pi*180
    angles[angles>180] = angles[angles>180]-360
    angles[angles<-180] = angles[angles<-180]+360
    
    lengths = np.absolute( gradmatrixx.flatten()+gradmatrixy.flatten()*1.j )
    nbins=np.r_[0:180:18]
    #axh.hist( angles[angles>0], weights=lengths[angles>0], bins=nbins, histtype='step' )
    #axh.hist( -angles[angles<0], weights=lengths[angles<0], bins=nbins, histtype='step' )
    
    radius = 2
    n=1
    A=int(np.pi*radius**2)
    nx, ny = n,n
    #nx, ny = np.array(matrix.shape)/divisions
    di, dj = np.array(matrix.shape)/n
    #print matrix.shape
    #print nx,ny
    xl = []
    ysig = []
    yamp = []
    nPix = len(ePix)
    meanE = np.sum(ePix)/nPix
    #di=0
    
    Xm, Ym = np.mgrid[:matrix.shape[0], :matrix.shape[1]]
    simulatedMatrix = np.zeros_like(matrix)
    
    simulatedMatrix2 = np.zeros_like(matrix)
    countMatrix = np.ones_like(matrix)
    ediffPix = np.zeros_like(ePix)


    tPix = np.dot( zip(xPix,yPix), v_[:,0] )
    lPix = np.dot( zip(xPix,yPix), v_[:,1] )
    tPix -= np.min(tPix)
    lPix -= np.min(lPix)
    meanT = np.average( tPix, weights=ePix )
    meanL = np.average( lPix, weights=ePix )
    stdT = std( tPix, weights=ePix )
    L = np.max(lPix)
    print 'props', meanT, meanL, stdT, L

    p12 = ( np.sum(ePix), meanT, stdT, 0)
    print 'initial guess', p12
    #func2 = lambda x,d,a,b,c: scipy.stats.norm.pdf( x[0]-a, 0, b + c*x[1] )*d
    func2 = lambda x, d, a, b, c: muonFit(x, 0, L, d, a, b, c)
    #fracp = -.2
    errfunc = lambda p_, x, y: ( func2( x, p_[0], p_[1], p_[2], p_[3] ) - y )*(2 + np.sign( func2( x, p_[0], p_[1], p_[2], p_[3] ) - y ) )
    pp3, success = scipy.optimize.leastsq(errfunc, p12[:], args=( np.vstack((tPix,lPix)), ePix) )
    print 'leastsq', pp3
    
    pp2 = scipy.optimize.curve_fit( func2, np.vstack((tPix,lPix)), ePix, p0=p12 )
    print 'first estimate', pp2[0]
    ediffPix = (ePix-func2(np.vstack((tPix,lPix)), *pp2[0]))**2
    mean_ediff, std_ediff = np.mean(ediffPix), np.std(ediffPix)
    print mean_ediff, std_ediff
    estCount = 1
    mask = np.abs(ediffPix-mean_ediff)<std_ediff
    looptPix, looplPix, loopePix = tPix[mask], lPix[mask], ePix[mask]
    while 1:
        estCount += 1
        pp22 = pp2[0][2]
        pp2 = scipy.optimize.curve_fit( func2, np.vstack((looptPix,looplPix)), loopePix, p0=pp2[0] )
        print estCount,'estimate', pp2[0]
        ediffPix = loopePix-func2(np.vstack((looptPix,looplPix)), *pp2[0])
        mean_ediff0, std_ediff0 = mean_ediff, std_ediff
        mean_ediff, std_ediff = np.mean(ediffPix), np.std(ediffPix)
        print mean_ediff, std_ediff
        err = np.abs(pp22-pp2[0][2])/np.abs(pp22)
        print 'err', err, len(loopePix), len(ePix)
        if err < 0.01: 
            break
        if estCount > 10:
            print 'did not converge'
            break
        #mask = np.abs(ediffPix-mean_ediff)<2*std_ediff
        mask = np.abs(ediffPix)<2*std_ediff
        #mask = ediffPix-mean_ediff<std_ediff
        #mask = loopePix < func2(np.vstack((looptPix,looplPix)), *pp2[0]) + 2*std_ediff
        looptPix, looplPix, loopePix = looptPix[mask], looplPix[mask], loopePix[mask]
    
    #simulatedMatrix2[ (xPix-np.min(xPix)).astype(int), (yPix-np.min(yPix)).astype(int) ] = func2(np.vstack((tPix,lPix)), *pp2[0])
    pedge = (L/2,L/2)
    pedge = scipy.optimize.curve_fit( lambda x,a,b: muonFit(x,a,b,*pp2[0]), np.vstack((looptPix,looplPix)), loopePix, p0=pedge )[0]
    pedge = (0,L)
    print 'edge', pedge
    #simulatedMatrix2[ (xPix-np.min(xPix)).astype(int), (yPix-np.min(yPix)).astype(int) ] = muonFit(np.vstack((tPix,lPix)), *[ pedge[0], pedge[1], pp2[0][0], pp2[0][1], pp2[0][2], pp2[0][3] ])
    simulatedMatrix2[ (xPix-np.min(xPix)).astype(int), (yPix-np.min(yPix)).astype(int) ] = muonFit(np.vstack((tPix,lPix)), *[ pedge[0], pedge[1], pp3[0], pp3[1], pp3[2], pp3[3] ])
    
    
    diffmatrix = (matrix-simulatedMatrix2)**2
    meandiff = np.mean(diffmatrix)
    stddiff = np.std(diffmatrix)
    maxdiff = np.max(diffmatrix)
    print meandiff, stddiff, maxdiff
    
    ax3 = fig.add_subplot(313)
    axh.imshow( (matrix-simulatedMatrix2).T, cmap=plt.cm.YlOrRd, origin='lower' )
    ax3.imshow( simulatedMatrix2.T, cmap=plt.cm.YlOrRd, origin='lower' )

    if 0:
        for i in range(di):
            for j in range(dj):
                mask = circularMask( xPix-np.min(xPix), yPix-np.min(yPix), i*nx, j*ny, radius )
                
                xPixij = xPix[mask].astype(float)
                yPixij = yPix[mask].astype(float)
                ePixij = ePix[mask]
                if len(ePixij) == 0: continue
                #print 'len', len(ePixij), np.sum(ePixij), np.sum(ePixij)/len(ePixij), meanE
                if np.sum(ePixij)/len(ePixij) < meanE: continue
                #print 'len', len(ePixij), A, .7*A
                #ax.add_patch(Circle( (i*nx, j*ny), radius, alpha=.25,edgecolor='black',facecolor='none'))
                mask2 = circularMask( X, Y, i*nx+1, j*ny+1, radius )
                mask3 = circularMask( Xm, Ym, i*nx+1, j*ny+1, radius )
                gmx = gradmatrixx[mask2]
                gmy = gradmatrixy[mask2]
                v0 = flow(0,0,gmx,gmy)
                v1 = flow(1,0,gmx,gmy)
                v2 = flow(2,0,gmx,gmy)
                v3 = flow(3,0,gmx,gmy)
                psi1 = np.angle(v1)
                psi2 = np.angle(v2)/2
                psi3 = np.angle(v3)/3
                v1_ = np.absolute(v1/v0)
                v2_ = np.absolute(v2/v0)
                v3_ = np.absolute(v3/v0)
                dirv2 = v2/np.absolute(v2)
                
                C = cov( np.vstack( (xPixij,yPixij) ), ePixij )
                w, v = np.linalg.eigh(C)
                angle0 = np.angle(v[0,1]+v[1,1]*1.j, deg=True)
                angle0rad = np.angle(v[0,1]+v[1,1]*1.j)
                angle1 = np.angle(v[0,0]+v[1,0]*1.j, deg=True)
                #ax.quiver( [i*nx], [j*ny], [v[0,1]], [v[1,1] ], units='dots', scale=.1, color='g', alpha=.5 )
                ax2 = fig2.add_subplot(di, dj, i*dj+(dj-j)+1)
                ax3 = fig3.add_subplot(di, dj, i*dj+(dj-j)+1)
                
                TV = lambda x: np.array([np.real(x),np.imag(x)])
                tPix = np.dot( zip(xPixij,yPixij), TV(dirv2) )
                t2Pix = np.dot( zip(xPixij,yPixij), v[:,0] )
                l2Pix = np.dot( zip(xPixij,yPixij), v[:,1] )
                dtPix = 1
                tPix -= np.min(tPix)
                t2Pix -= np.min(t2Pix)
                l2Pix -= np.min(l2Pix)
                N=30
                dT = (np.max(tPix)-np.min(tPix))/N
                #ax3.errorbar( tPix, ePixij, xerr=0, fmt='g.', markersize = .1 )
                ax3.errorbar( t2Pix, ePixij, xerr=0, fmt='r.', markersize = .1 )
                p1 = (np.average( t2Pix, weights=ePixij ), std( t2Pix, ePixij ), np.sum(ePixij) )
                
                pp = None
                try:
                #if 1:
                    func = lambda x,a,b,c: scipy.stats.norm.pdf(x,a,b)*c
                    pp = scipy.optimize.curve_fit( func, t2Pix, ePixij, p0=p1 )
                    ax3.plot( np.array(range(N))*dT, func( np.array(range(N))*dT, *pp[0]), 'b-' )
                    ax3.tick_params(labelbottom='off',labelleft='off')
                    
                    tcenter = pp[0][0]-np.mean(t2Pix)
                    #print radius, pp[0][0]-np.mean(t2Pix)
                    ang = np.angle( v[0,0]+v[1,0]*1.j )
                    simulatedMatrix[ (xPixij-np.min(xPix)).astype(int), (yPixij-np.min(yPix)).astype(int) ] = func(t2Pix, *pp[0])
                    countMatrix[ (xPixij-np.min(xPix)).astype(int), (yPixij-np.min(yPix)).astype(int) ] += 1
                    
                    
                    #if abs( tcenter ) < .75*radius:
                        #ax.plot( 
                            #[i*nx + radius*np.cos(ang+np.pi/2+np.pi) + tcenter*np.cos(ang),
                            #i*nx + radius*np.cos(ang+np.pi/2) + tcenter*np.cos(ang)], 
                            #[j*ny + radius*np.sin(ang+np.pi/2+np.pi) + tcenter*np.sin(ang), 
                            #j*ny + radius*np.sin(ang+np.pi/2) + tcenter*np.sin(ang)],  'b-')
                except:
                    #print sys.exc_info()[0]
                    pass

                
                #ax.quiver( [i*nx,i*nx], [j*ny,j*ny], [v2_*np.cos(psi2),v2_*np.cos(psi2+np.pi)], [v2_*np.sin(psi2),v2_*np.sin(psi2+np.pi)], units='dots', scale=.05, color='r' )
                #ax.quiver( 
                    #[i*nx,i*nx,i*nx], 
                    #[j*ny,j*ny,j*ny], 
                    #[v3_*np.cos(psi3),v3_*np.cos(psi3+2*np.pi/3),v3_*np.cos(psi3+4*np.pi/3)], 
                    #[v3_*np.sin(psi3),v3_*np.sin(psi3+2*np.pi/3),v3_*np.sin(psi3+4*np.pi/3)], units='dots', scale=.05, color='g' )
                
                angles = np.angle( gmx.flatten()+gmy.flatten()*1.j, deg=True )
                angles = angles - psi2/np.pi*180
                maxAngle = 180
                angles[angles>180] = angles[angles>180]-360
                angles[angles<-180] = angles[angles<-180]+360
                angles[angles>90] = 180-angles[angles>90]
                angles[angles<-90] = -180-angles[angles<-90]
                lengths = np.absolute( gmx.flatten()+gmy.flatten()*1.j )
                nbins=np.r_[0:91:90./5]
                #for a,l in zip(angles[angles>0],lengths[angles>0]):
                    #axh.arrow( a, 0, 0, l ) 
                ##for a,l in zip(angles[angles<0],lengths[angles<0]):
                    ##axh.arrow( a+np.pi/2, 0, 0, l ) 
                ax2.hist( angles[angles>0], weights=lengths[angles>0], bins=nbins, histtype='step' )
                ax2.hist( -angles[angles<0], weights=lengths[angles<0], bins=nbins, histtype='step' )
                ax2.tick_params(labelbottom='off',labelleft='off') 
                if not pp is None:
                    if np.abs(tcenter)<radius:
                        xl += [ np.dot( (i*nx,j*ny), v_[:,1] ) ]
                        ysig += [ pp[0][1] ]
                        yamp += [ pp[0][2] ]

    #axh.imshow( (simulatedMatrix/countMatrix).T, cmap=plt.cm.YlOrRd, origin='lower' )
    #if len(ysig)>0:
        #ax = fig.add_subplot(313)
        #ax.plot((np.array(xl)-np.min(xl))/(np.max(xl)-np.min(xl)), np.array(ysig)/np.max(ysig), 'b.')
        #ax.plot((np.array(xl)-np.min(xl))/(np.max(xl)-np.min(xl)), np.array(yamp)/np.max(yamp), 'r.')
    
    fig.subplots_adjust()
    fig.savefig('viewer.pdf')
    fig2.subplots_adjust()
    fig2.savefig('viewer2.pdf')
    fig3.subplots_adjust()
    fig3.savefig('viewer3.pdf')
    return
    
    
def plotHist(fname, expr, sel, Max = None, Min = None, expr2 = None, sel2 = None ):
    baseSel = 'flag==0' + ('&&%s<%s'%(expr,Max) if Max else '') + ('&&%s>%s'%(expr,Min) if Min else '')
    dataAll = root_numpy.root2array(fname, treename = 'hitSumm', branches = expr, selection = baseSel )
    dataSel = root_numpy.root2array(fname, treename = 'hitSumm', branches = expr, selection = baseSel+('&&'+sel if not sel == "" else "") )
    dataSel2 = None
    data2Sel = None
    data2Sel2 = None
    if not sel2 is None:
        dataSel2 = root_numpy.root2array(fname, treename = 'hitSumm', branches = expr, selection = baseSel+('&&'+sel2 if not sel2 == "" else "") )
    fig = plt.figure(figsize=(10,10))
    fig.subplots_adjust(hspace=0, wspace=0)
    
    if not expr2 is None:
        ax = fig.add_subplot(221)
        dataAll2 = root_numpy.root2array(fname, treename = 'hitSumm', branches = expr2, selection = baseSel )
        data2Sel = root_numpy.root2array(fname, treename = 'hitSumm', branches = expr2, selection = baseSel+('&&'+sel if not sel == "" else "") )
        if not sel2 is None: data2Sel2 = root_numpy.root2array(fname, treename = 'hitSumm', branches = expr2, selection = baseSel+('&&'+sel2 if not sel2 == "" else "") )
        
    else:
        ax = fig.add_subplot(111)
    ax.set_title( expr )
    print dataSel.shape
    _, bins, _ = ax.hist( dataSel, bins = int(np.sqrt(len(dataSel))), histtype='step', label=sel, log=True )
    ax.hist( dataAll, bins, histtype='step', label='all', log=True )
    if not sel2 is None: ax.hist( dataSel2, bins, histtype='step', label=sel2, log=True )
    ax.legend( fancybox=True, framealpha=0.1)
    ax.grid(True)
    if not expr2 is None:
        ax2 = fig.add_subplot(224)
        ax2.set_xlabel( expr2 )
        _, bins, _ = ax2.hist( data2Sel, bins = int(np.sqrt(len(data2Sel))), histtype='step', label=sel, log=True, orientation = "horizontal" )
        ax2.hist( dataAll2, bins, histtype='step', label='all', log=True, orientation = "horizontal" )
        if not sel2 is None: ax2.hist( data2Sel2, bins, histtype='step', label=sel2, log=True, orientation = "horizontal" )
        ax2.legend( fancybox=True, framealpha=0.1)
        ax2.grid(True)
        
        ax3 = fig.add_subplot(223, sharey=ax2, sharex=ax)
        print data2Sel.shape
        ax3.plot( dataSel, data2Sel, '.', ms=.1 )
        if not sel2 is None: ax3.plot( dataSel2, data2Sel2, '.', ms=.1 )
        ax3.grid(True)
        
    fig.savefig('viewer.pdf')
    return

if __name__ == "__main__":
    print "viewer by Philipe Mota (philipe.mota@gmail.com)"
    args = sys.argv
    if '-h' in args:
        print usagestr
        exit(0)
    f       = None
    p = False
    i = 0
    expr    = None
    expr2    = None
    logy     = False
    sel = ""
    sel2 = None
    lim = None
    Max = None
    Min = None
    imin = None
    imax = None
    scan = False
    l = False
    if '-f' in args:
        try: 
            f = args[ args.index( '-f' ) + 1 ]
        except:
            print "no input file given."
            exit(0)
        print 'input file', f
    if '-l' in args:
        l = True
    if '-p' in args:
        p = True
    if '-scan' in args:
        scan = True
    if '-i' in args:
        try: 
            i = int(args[ args.index( '-i' ) + 1 ])
        except:
            print "no expression given."
        print 'index', i
    if '-imin' in args:
        try: imin = int(args[ args.index( '-imin' ) + 1 ])
        except: print "no expression given."
        print 'indexmin', imin
    if '-imax' in args:
        try: imax = int(args[ args.index( '-imax' ) + 1 ])
        except: print "no expression given."
        print 'indexmax', imax
    if '-x' in args:
        try: 
            expr = args[ args.index( '-x' ) + 1 ]
        except:
            print "no expression given."
        print 'xAxis is', expr
    if '-x2' in args:
        try: 
            expr2 = args[ args.index( '-x2' ) + 1 ]
        except:
            print "no expression given."
        print 'xAxis is', expr
    if '-max' in args:
        try: 
            Max = float(args[ args.index( '-max' ) + 1 ])
        except:
            print "no expression given."
        print 'max limit', Max
    if '-min' in args:
        try: 
            Min = float(args[ args.index( '-min' ) + 1 ])
        except:
            print "no expression given."
        print 'min limit', Min
    if '-s' in args:
        try: 
            sel = args[ args.index( '-s' ) + 1 ]
        except:
            print "no selection given."
        print 'selection is', sel
    if '-s2' in args:
        try: 
            sel2 = args[ args.index( '-s2' ) + 1 ]
        except:
            print "no selection given."
        print 'selection is', sel
    if '-logy' in args:
        logy = True
        
    if l:
        print 'list of branches in hitSumm'
        print root_numpy.list_branches(f, treename='hitSumm')
        print
    if f:
        if scan:
            scanTracks( f, imin, imax )
            exit(0)
        if p:
            plotTrack(f, index=i)
            exit(0)
        plotHist(f, expr, sel, Max = Max, Min = Min, expr2 = expr2, sel2 = sel2)
    exit(0)
