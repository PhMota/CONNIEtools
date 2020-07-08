from __future__ import print_function
import time
import datetime
from termcolor import colored

class Timer:
    def __init__(self, msg='elapsed'):
        self.msg = msg
        
    def __enter__(self):
        self.start = datetime.datetime.now()
        return self
    
    def __del__(self):
        s = [ self.msg ]
        s += [ colored(self.time(),'green') ]
        print( ' '.join(s) )
    
    def __exit__(self, type, value, traceback):
        s = [ self.msg ]
        if value:
            s += [ colored(self.time(),'red') ]
            s += [ '%s'%type ]
            s += [ '%s'%value ]
        else:
            s += [ colored(self.time(),'green') ]
        print( ' '.join(s) )
    
    def seconds(self, N=1):
        return (datetime.datetime.now() - self.start).total_seconds()*N
    
    def time(self, N=1):
        t = self.seconds(N)
        return self.make_str_ms(t)

    def make_str_ms(self, t):
        if t > 60*60:
            s = '%sh%sm%ss%sms' %( int(t)/60/60, int(t/60)%60, int(t)%60, int(t*1e3)%1000 )
        else: 
            if t > 60: s = '%sm%ss%sms' %( int(t)/60, int(t)%60, int(t*1e3)%1000 )
            else: s = '%ss%sms' % (int(t), int(t*1e3)%1000 )
        return s

    def make_str_s(self, t):
        if t > 60*60:
            s = '%sh%sm%ss' %( int(t)/60/60, int(t/60)%60, int(t)%60 )
        else: 
            if t > 60: s = '%sm%ss' %( int(t)/60, int(t)%60 )
            else: s = '%ss' % (int(t) )
        return s
        
    def check(self, wait_secs=None, total_loop_number=None):
        if not hasattr(self,'wait_secs'): self.wait_secs = wait_secs
        if not hasattr(self, 'total_loop_number'): self.total_loop_number = total_loop_number
        try:
            self.loop_number +=1
        except AttributeError:
            print('start loop timer:', 'next report in {}'.format( self.make_str_s(self.wait_secs) ), 'total', self.total_loop_number )
            self.loop_number = 0
            self.count = 1
            return
        if self.seconds()/self.wait_secs > self.count:
            secs_per_loop = self.seconds()/self.loop_number
            remaining_loops = self.total_loop_number - self.loop_number
            eta_secs = secs_per_loop * remaining_loops
            if eta_secs < self.wait_secs and eta_secs > 30:
                self.wait_secs /= 2
                self.count *= 2
            print('eta', colored( self.make_str_s( eta_secs ), 'yellow'), 'done', self.loop_number, '[{:.2g}]'.format(float(self.loop_number)/self.total_loop_number), 'next report in {}'.format( self.make_str_s(self.wait_secs) ) )
            self.count += 1
            if eta_secs > 3*self.wait_secs: 
                self.wait_secs *= 2
                self.count /= 2
        return
    
    def _wait_secs(self, wait_secs):
        self.wait_secs = wait_secs
        
    def _total_loop_number(self, tln):
        self.total_loop_number = tln
    
    def __next__(self):
        return self.check()
    next = __next__
    
    def __iter__(self):
        return self

class timer:
    def __init__(self, name = ''):
        self.name = name
        self.start = time.clock()
        self.text = ''
    def factor( self, text = '', f = None ):
        self.text = text
        self.f = f
    def __del__(self):
        print( self.__call__() )
        if not self.text == '': print( '\t%s: %.3gs'%(self.text, self.f*elapsed) )
    def __call__(self):
        elapsed = time.clock() - self.start
        return 'timer[%s] '%self.name + str(datetime.timedelta(seconds=elapsed))

class LoopTimer:
    def __init__(self, n=1, name=''):
        self.name = name
        self.start = time.time()
        self.n = n
        print( self.name, 'start', time.strftime("%Hh%Mm%Ss", time.localtime(self.start) ) )
        return
    
    def __del__(self):
        print( self.name, 'done in', time.strftime("%Hh%Mm%Ss", time.gmtime(time.time() - self.start) ) )
        return
    
    def eta( self, i ):
        if i==0: return
        ave = (time.time() - self.start)/i
        print( self.name, 'eta', time.strftime("%Hh%Mm%Ss", time.gmtime((self.n-i)*ave ) ) )
        return
    
    def end( self, i ):
        if i==0: return
        ave = (time.time() - self.start)/i
        print( self.name, 'end', time.strftime("%Hh%Mm%Ss", time.localtime(time.time() + (self.n-i)*ave ) ) )
        return
        
