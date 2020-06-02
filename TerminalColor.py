
def cmd(x):
    return '\033[%dm' % x

Q = cmd(0)
B = cmd(1)
I = cmd(3)
U = cmd(4)
r = cmd(91)
g = cmd(32)
y = cmd(33)
b = cmd(34)

colorscaleYR = staticmethod( lambda x, min_, max_: '\033[48;5;%dm \033[0m' % (1+6*int( (( 196/6 - 226/6 )*( x - min_)/( max_ - min_ ) + 226/6 ))) )
colorscaleGray = staticmethod( lambda x, min_, max_: ('\033[48;5;%dm \033[0m' % (int( (( 255 - 232 )*( x - min_)/( max_ - min_ ) + 232 ))) )*2 )

def text( x, color=None, mode=None ):
    colorf = lambda x: x
    if color == 'red' or color == 'r':
        colorf = lambda x: r + x + Q
    elif color == 'blue' or color == 'b':
        colorf = lambda x: b + x + Q
    elif color == 'green' or color == 'g':
        colorf = lambda x: g + x + Q
    elif color == 'yellow' or color == 'y':
        colorf = lambda x: y + x + Q
    
    modef = lambda x: x
    if mode == 'bold' or mode == 'B':
        modef = lambda x: B + x + Q
    elif mode == 'italic' or mode == 'I':
        modef = lambda x: I + x + Q
    elif mode == 'underline' or mode == 'U':
        modef = lambda x: U + x + Q
    
    return modef( colorf(x) )
