from __future__ import print_function
import subprocess
import os
from fstring import F
from termcolor import colored
from progress import progressbar

preamble = r'''
\documentclass{beamer}
%%
%% Choose how your presentation looks.
%%
%% For more themes, color themes and font themes, see:
%% http://deic.uab.es/~iblanes/beamer_gallery/index_by_theme.html
%%
\mode<presentation>
{
  \usetheme{default}      %% or try Darmstadt, Madrid, Warsaw, ...
  \usecolortheme{default} %% or try albatross, beaver, crane, ...
  \usefonttheme{default}  %% or try serif, structurebold, ...
  \setbeamertemplate{navigation symbols}{}
  \setbeamertemplate{caption}[numbered]
  \setbeamertemplate{footline}[frame number]
}

\usepackage[english]{babel}
\usepackage[utf8x]{inputenc}
\usepackage{grffile}
\usepackage{hyperref}
\usepackage{listings}
\usepackage{xcolor}
\usepackage{nath}
\delimgrowth=1
\setlength{\mathindent}{-1in}

\usefonttheme{professionalfonts}

\lstset{
language=Python,
breakatwhitespace=true,
breaklines=true,
%% basicstyle=\fontsize{4}{5}\sffamily,
basicstyle=\fontsize{8}{10}\sffamily,
keywordstyle=\bf\color{blue},
showtabs=true,
tabsize=2,
%% digitstyle=\color{red},
%% otherstyle=\color{red},
numberstyle=\color{red},
columns=fullflexible,
%% commentstyle=\sf\tiny\color{gray},
%% stringstyle=\color{red},
showstringspaces=false,
morekeywords={as},
%% emph={(,)},
%% emphstyle=\color{blue},
numbers=left
%% identifierstyle=\color{blue}
}
\lstset{literate={-}{-}1}
\title[%s]{%s}
\author{Philipe Mota}
\institute{CBPF}
\date{\today}

\begin{document}

\begin{frame}
  \titlepage
\end{frame}
'''

B = '\n'

#class Frame:

class Beamer:
    def writelines( self, mode = 'a' ):
        tab = '\t'
        endl = '\n'
        level = 0
        with open( '%s.tex' % (self.folder+'/'+self.basename), mode ) as file:
            for element in progressbar( self.elements, msg='document' ):
                if callable( element ):
                    element = element()
                if not hasattr( element, '__iter__' ):
                    element = [ element ]
                for e in element:
                    e = str(e)
                    n = 1
                    if r'\end' in e or r'\section' in e or r'\part' in e:
                        n = 2
                    if r'\end' in e:
                        level -= 1
                    try:
                        file.write( tab*level + e + endl*n )
                    except TypeError:
                        print( 'failed to write', e, type(e) )
                        raise TypeError
                    if r'\begin' in e:
                        level += 1

    def __init__(self, folder='calculations', title='simulation tools'):
        self.folder = folder
        if not os.path.exists(self.folder):
            os.mkdirs(self.folder)
        self.basename = self.folder.split('/')[-1]
        self.elements = []
        self.elements.append( preamble %( title, title ) )

    def __enter__(self):
        return self

    def __exit__(self, type, value, traceback):
        self.elements.append(r'\end{document}')
        self.writelines( mode = 'w' )
        self.pdflatex()

    def pdflatex(self):
        self.subprocess_cmd( 'cd {0}; pwd; pdflatex {1}.tex >/dev/null; cd ..'.format(self.folder, self.basename) )

    def subprocess_cmd(self, command):
        subprocess.call( [command], shell=True )

    def par( self, s ):
        self.elements.append( s )

    def __iadd__( self, elements ):
        # print( '__iadd__', elements )
        if not hasattr(elements, '__iter__'):
            elements = [elements]
        self.elements.extend( elements )
        return self

    equation_evironment = 'equation'

    def section( self, text ):
        self.elements.append( [ r'\section{%s}' % text] )
        # self.writelines( [ r'\section{%s}' % text] )

    def environment( self, s, env = 'equation' ):
        # self.writelines( [r'\begin{%s}' % env] )
        # self.writelines( s )
        # self.writelines( [r'\end{%s}' % env] )
        self.elements.append( [r'\begin{%s}' % env] )
        self.elements.append( s )
        self.elements.append( [r'\end{%s}' % env] )

    def eq( self, lhs, rhs = None, mathmode = 'equation' ):
        #print 'eq', lhs
        r = [r'\begin{%s}'%mathmode, lhs]
        if not rhs is None:
            if type(rhs) is list or type(rhs) is tuple:
                r += ['=', r'\\ '.join(rhs)]
            else:
                r += ['=', rhs]
        r += [r'\end{%s}'%mathmode]
        # self.writelines( r )
        self.elements.append(r)

    def math( self, lhs, rhs = None, label = None ):
        #print 'eq', label
        r = [r'\begin{equation}', latex(lhs)]
        if not rhs is None:
            r += ['=', latex(rhs)]
        if not label is None:
            r += [r'\label{%s}'%label]
        r += [r'\end{equation}']
        # self.writelines( r )
        self.elements.append(r)

    def matheval( self, s ):
        #print 'eq', label
        # self.writelines( [r'\begin{equation}', latex(eval(m)), r'\label{%s}'%m, r'\end{equation}'] )
        self.elements.append( [r'\begin{equation}', latex(eval(m)), r'\label{%s}'%m, r'\end{equation}'] )

    def center(self, *lines ):
        s = r'\begin{center}' + B
        for line in lines:
            s += line
        s += r'\end{center}' + B
        return s

    def set_func( self, func ):
        self.func = func

    def check_file( self, fname ):
        if not os.path.exists(fname):
            print('!!! file not found', fname)
            print('running function' )
            self.func()
            if not os.path.exists(fname):
                print('!!! file not generated', fname)
                exit()
        print('file found', fname)
        print('skipping')

    def table( self, file, fontsize=12, spacing=15, func=None, divide=1 ):
        s = ''
        fname = '%s/%s' % (self.fname,file)
        self.check_file(fname)

        for d in range(divide):
            s += r'{{\fontsize{{{}}}{{{}}}\selectfont'.format(fontsize,spacing)+B
            s += r'{\setlength{\tabcolsep}{0.2em}'+B
            for line in open(fname, 'r'):
                args = line.strip('\n').split(',')
                args = args[ d*len(args)/divide:(d+1)*len(args)/divide ]
                if '#' in line:
                    args[0] = args[0].strip('#')
                    s += r'\begin{{tabular}}{{ {} }}'.format('r'*len(args)) +B
                s += ' & '.join( args ) + r'\\' + B
            s += r'\end{tabular}'+B
            s += r'}'+B
            s += r'}'+B
        return s

    def tabular( self, matrix, align=None, s='' ):
        if align in ['c','r','l']:
            align = [align]*len(matrix[0])
        s += r'{\setlength{\tabcolsep}{0em}'+B
        s += r'\begin{center}' + B
        s += r'\begin{{tabular}}{{ {} }}'.format( ''.join(align) ) + B
        s += '\\\\ \n'.join( [ ' & '.join( line ) for line in matrix ] )
        s += r'\end{tabular}' + B
        s += r'\end{center}' + B
        s += r'}' + B
        #print( s )
        return s

def textcolor(c,s):
    return r'\textcolor{%s}{ %s }' % (c,s)

def part(name):
    return [ F(r'\part{ {{name}} }'), F(r'\frame{\partpage}') ]

def figure( path, width=None, height=None, scale=None, center=True ):
    elements = []
    fname = path
    if center:
        elements.append( r'\begin{center}' )
    if width is not None:
        elements.append(
            r'\includegraphics[width={width}\columnwidth]{{{fname}}}'.format(width=width, fname=fname)
        )
    elif height is not None:
        elements.append(
            r'\includegraphics[height={height}\paperheight]{{{fname}}}'.format(height=height, fname=fname)
        )
    elif scale is not None:
        elements.append(
            r'\includegraphics[scale={scale}]{{{fname}}}'.format(scale=scale, fname=fname)
        )
    else:
        raise Exception('missing heigh or width')
    if center:
        elements.append( r'\end{center}' )
    return elements

def code( c=None, language=None ):
    return [
        r'\begin{lstlisting}[language=%s]' % (language),
        c,
        r'\end{lstlisting}'
    ]

def inlinecode( c, language='Python' ):
    return '\lstinline[language=%s]{%s}' % (language,c)

def frame( title, *args ):
    elements = []
    elements.append( r'\begin{frame}[fragile]{%s}' % title )
    for arg in args:
        if not hasattr(arg, '__iter__'): arg = [arg]
        elements.extend( arg )
    elements.append( r'\end{frame}' )
    return elements

def column( *args, **kwargs ):
    elements = []
    widths = kwargs['widths']
    elements.append( r'\begin{columns}' )
    for arg, width in zip(args, widths):
        elements.append(
            r'\column{%s\paperwidth}'%( width if width>0 else -sum( widths ) )
        )
        if not hasattr(arg, '__iter__'):
            arg = [arg]
        elements.extend( arg )
    elements.append( r'\end{columns}' )
    return elements

def itemize( *items ):
    elements = []
    elements.append(r'\begin{itemize}')
    for item in items:
        elements.append('\item ' + item)
    elements.append(r'\end{itemize}')
    return elements

def openBeamer( fname, title ): return Beamer( fname, title )

#def math( expr ):
    #return '\n'.join( ['$$', expr, '$$'] )

#def __invert__(a):
    #return 'test'

class MathGen:
    def __or__(self, a):
        return '\n'.join( ['$$', str(a), '$$'] )

#Math = MathGen()
def Math( a ):
    if type(a) is tuple or type(a) is list:
        return '\n'.join( ['\n$$', ',\quad '.join(map(str,a)), '.$$\n'] )
    return '\n'.join( ['\n$$', str(a), '.$$\n'] )

def delim( a ):
    if '(' in str(a):
        if '[' in str(a):
            return r'\{{ {} \}}'.format( str(a) )
        return Expr( r'[{}]'.format( str(a) ) )
    return Expr( F(r'({{str(a)}})') )

class Expr:
    def __init__(self, s):
        self.s = s

    def __str__(self):
        return str(self.s)

    def __repr__(self):
        return str(self.s)

    def __neg__(self):
        return Expr(r'-{}'.format(self.s))

    def __add__(self, a):
        return ExprAdd(r'{}+{}'.format(self.s, str(a)))

    def __radd__(self, a):
        return ExprAdd(r'{}+{}'.format( str(a), self.s ))

    def __sub__(self, a):
        if isinstance(a, ExprAdd):
            a = delim(a)
        return ExprAdd(r'{}-{}'.format(self.s, str(a)))

    def __rsub__(self, a):
        return ExprAdd(r'{}-{}'.format( str(a), self.s ))

    def __div__( self, a ):
        return frac(self.s, a)

    def __rdiv__( self, a ):
        return frac(a, self.s )

    def __mul__( self, a ):
        if isinstance( a, ExprAdd ):
            return a.__rmul__(self.s)
        return ExprMul(r'{}{}'.format(self.s, str(a) ))

    def __rmul__( self, a ):
        return ExprMul(r'{}{}'.format( str(a), self.s ))

    def __or__( self, a ):
        return Expr(r'{}|{}'.format(self.s, str(a) ))

    def __call__( self, *a ):
        p = ', '.join( map(str,a) )
        p = delim(p)
        return Expr( r'{}{}'.format( self.s, p ) )

    def __eq__( self, a ):
        if type(a) is tuple or type(a) is list:
            return Expr(r'{} \wall = {} \return'.format( self.s, r' \\ ='.join(map(str,a)) ) )
        return Expr(r'{} = {}'.format( self.s, str(a) ) )

    def __req__( self, a ):
        return Expr(r'{} = {}'.format( str(a), self.s ) )

    def __getitem__( self, a ):
        if type(a) is tuple:
            return Expr(r'{}_{{{}}}^{{{}}}'.format( self.s, str(a[0]), str(a[1]) ) )
        return Expr(r'{}_{{{}}}'.format( self.s, str(a) ) )

    def __pow__( self, a ):
        return Expr(r'{}^{{ {} }}'.format( self.s, str(a) ) )

    def __rpow__( self, a ):
        return Expr(r'{}^{{ {} }}'.format( str(a), self.s ) )

    def __sqrt__(self):
        return Expr(r'\sqrt{{ {} }}'.format( self.s ) )

class ExprGen:
    def __or__( self, a ):
        return Expr(a)

class ExprAdd(Expr):
    def __init__(self, s):
        self. s = s

    def __pow__( self, a ):
        return Expr(r'{}^{}'.format( delim(self.s), str(a) ) )

    def __mul__( self, a ):
        return ExprMul(r'{}{}'.format( delim(self.s), str(a) ) )

    def __rmul__( self, a ):
        return ExprMul(r'{}{}'.format( str(a), delim(self.s) ) )

    def __neg__(self):
        return Expr(r'-{}'.format(delim(self.s)))

    def __rsub__(self, a):
        return Expr(r'{}-{}'.format(str(a),delim(self.s)) )


class ExprMul(Expr):
    def __init__(self, s):
        self. s = s
    def __pow__( self, a ):
        return Expr(r'{}^{{ {} }}'.format( delim(self.s), str(a) ) )


mu = Expr(r'\mu ')
nu = Expr(r'\nu ')
delta = Expr(r'\delta ')
sigma = Expr(r'\sigma ')
alpha = Expr(r'\alpha ')
pi = Expr(r'\pi ')
xi = Expr(r'\xi ')

Gamma = Expr(r'\Gamma ')
Delta = Expr(r'\Delta ')
Prod = Expr(r'\prod ')
Sum = Expr(r'\sum ')
exp = Expr(r'{\rm exp}')
e = Expr(r'{\rm e}')
log = Expr(r'{\rm log}')
erf = Expr(r'{\rm erf}')

Int = Expr(r'\int')
inf = Expr(r'\infty')
oo = Expr(r'\infty')

E = ExprGen()

def frac(a, b):
    return Expr( r'\frac{{ {} }}{{ {} }}'.format( str(a), str(b)) )

def bf(a):
    return Expr( r'{{\mathbf {} }}'.format(str(a)) )

def ave(a):
    return Expr( r'\langle {} \rangle '.format(str(a)) )

def bar(a):
    return Expr( r'\bar{{ {} }}'.format(str(a)) )

def hat(a):
    return Expr( r'\hat{{ {} }}'.format(str(a)) )

def sqrt(a):
    return Expr(r'\sqrt{{ {} }}'.format( str(a) ) )

def d(a):
    return Expr(r'{{\rm d}} {}\, '.format( str(a) ) )

def factorial(a):
    if isinstance(a, ExprAdd) or isinstance(a, ExprMul):
        Expr(r'({})!'.format( str(a) ) )
    return Expr(r'{}!'.format( str(a) ) )


import hashlib

def my_hash(entry):
    return int(hashlib.md5(str(entry).encode('utf-8')).hexdigest(), 16)

def makefigure( code, filename, height=1., folder='', nocode=False ):
    split_filename = str(filename).split('.')
    print( 'hash', str(my_hash(code)) )
    filename = '.'.join( split_filename[:-1] + ['_hash'+str(my_hash(code)), split_filename[-1]] )
    fullpath = F('{{folder}}/{{filename}}').str()
    print( 'fullpath', fullpath )

    code += F(' --output "{{fullpath}}"').str()
    code_print = code.replace('\n', r'\n')
    elements = []
    if not nocode:
        elements.extend( code( code_print, language='Bash' ) )
    if not os.path.exists(fullpath):
        # print()
        # print( colored(code_print, 'green' ))
        # print()
        # subprocess.call( ['%s' % code], shell=True )
        code_split = code.split(' ')
        cmd = ' '.join(code_split[:2])
        # print( 'code_split', cmd )
        args = ' '.join(code_split[2:])
        # print( 'code_split', args )
        subprocess.call( F('{{cmd}} {{args}}').str(), shell=True )
    if not os.path.exists(fullpath):
        print()
        print( '!!!not found', fullpath )
        exit(0)
    elements.extend( figure( filename, height=height, center=True ) )
    return elements

def requiredFile( code, file, doc ):
    a = doc.code( code, language='Bash' )
    if not os.path.exists(folder+'/'+file):
        print()
        print(code)
        print()
        subprocess.call( [code], shell=True )
    #b = doc.figure( file, height=height, center=True )
    return a

def makecmd( **options ):
    cmd = options.pop('prog')
    print( options )
