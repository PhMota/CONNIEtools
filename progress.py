from __future__ import print_function
import sys
import time
from termcolor import colored

class progressbar:
    active = []
    def __init__(self, iter, msg='', length=50):
        self.iter = iter
        self.current = 0
        self.length = length
        self.msg = msg
        progressbar.active.append(self)
    def iter(self):
        return self
    __iter__ = iter
    def next(self):
        self.current += 1
        for i, active in enumerate(progressbar.active[::-1]):
            complete = 100*(active.current-1)/len(active.iter)
            print( colored('|' + '#'*(complete*active.length/100) + ' '*(active.length -complete*active.length/100) + '| %3s%%' % complete + ' %s' %active.msg, 'green') )

        for active in progressbar.active:
            sys.stdout.write("\033[F") # Cursor up one line
            sys.stdout.write("\033[K") # Clear to the end of line

        try:
            return self.iter[self.current-1]
        except IndexError:
            progressbar.active.pop()
            raise StopIteration


if __name__ == '__main__':
    for i in progressbar(range(5)):
        print(i)
        print(i+1)
        print(i**2)
        time.sleep(1)

        for j in progressbar(range(10)):
            time.sleep(1)
            print(j)
            print(j+1)
            print(j**2)
        for j in progressbar(range(10)):
            time.sleep(1)
            print(i*j)
            print(i+1)
            print(i**2)

    for i in progressbar(range(5)):
        sys.stdout.write("\033[K") # Clear to the end of line
        print(i)
        print(i+1)
        print(i**2)
        print("Loading" + "." * i)
        sys.stdout.write("\033[F") # Cursor up one line
        time.sleep(1)

    for i in range(5):
        print("Loading" + "." * i)
        sys.stdout.write("\033[F") # Cursor up one line
        time.sleep(1)
