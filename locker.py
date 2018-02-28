import os
import time

class Locker:
    def __init__( self, user, fname ):
        self.user = user
        self.fname = fname
        if os.path.exists( self._makename_() ):
            print 'sleeping: locker found', self._makename_()
        while os.path.exists( self._makename_() ):
            time.sleep(1)
        print 'lock', self._makename_()
        open( self._makename_(), 'w' ).close()
        return
        
    def _makename_(self):
        return "%s.%s.lock"%(self.fname, self.user)

    def free( self ):
        if os.path.exists( self._makename_() ):
            os.remove( self._makename_() )
        print 'locker', self._makename_(), 'removed'
        return
        
    def __del__( self ):
        self.free()
        return