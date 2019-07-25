import time
import datetime

class timer:
    def __init__(self, name = ''):
        self.name = name
        self.start = time.clock()
        self.text = ''
    def factor( self, text = '', f = None ):
        self.text = text
        self.f = f
    def __del__(self):
        print self.__call__()
        if not self.text == '': print '\t%s: %.3gs'%(self.text, self.f*elapsed)
    def __call__(self):
        elapsed = time.clock() - self.start
        return 'timer[%s] '%self.name + str(datetime.timedelta(seconds=elapsed))

class LoopTimer:
    def __init__(self, n=1, name=''):
        self.name = name
        self.start = time.time()
        self.n = n
        print self.name, 'start', time.strftime("%Hh%Mm%Ss", time.localtime(self.start) )
        return
    def __del__(self):
        print self.name, 'done in', time.strftime("%Hh%Mm%Ss", time.gmtime(time.time() - self.start) )
        return
    def eta( self, i ):
        if i==0: return
        ave = (time.time() - self.start)/i
        print self.name, 'eta', time.strftime("%Hh%Mm%Ss", time.gmtime((self.n-i)*ave ) )
        return
    def end( self, i ):
        if i==0: return
        ave = (time.time() - self.start)/i
        print self.name, 'end', time.strftime("%Hh%Mm%Ss", time.localtime(time.time() + (self.n-i)*ave ) )
        return
        
