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
