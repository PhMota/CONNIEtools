from __future__ import print_function
import os
from termcolor import colored

rows, columns = os.popen('stty size', 'r').read().split()
def print_var( var, vars_, end = '\n', line_char='', start_pos=0 ):
    if not type(var) is str:
        start_pos = 0
        for ivar in var:
            start_pos = print_var( ivar, vars_, end=' ', start_pos=start_pos, line_char=line_char ) + 1
        print()
    else:
        val = vars_[var]
        if type(val) is float: val = '{:.4f}'.format(val)
        text = '%s %s' % ( colored(var, 'green' ), val )
        
        #print(start_pos)
        if start_pos == 0:
            text = line_char + text
        elif start_pos + len(text) >= int(columns):
            start_pos = 0
            print()
            text = line_char + text

        print( text, end=end )
        return start_pos + len(text)
