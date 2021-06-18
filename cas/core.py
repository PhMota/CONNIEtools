import IPython.display

def display(s, eval=False):
    if eval == True:
        IPython.display.display( IPython.display.Math( fr"{s} = {s.eval()}" ) )
    else:
        IPython.display.display( IPython.display.Math( fr"{s}" ) )

def delim(s):
    return fr"\left({s}\right)"

def Symbol(name):
    return SymbolClass(name)

def Add(*args):
    return AddClass(*args)

def Any(name):
    return AnyClass(name)

def Delta(arg):
    return DeltaClass(arg)

class SumMaker(object):
    def __getitem__(self, lims):
        def factory(arg):
            return SumClass(lims, arg)
        return factory

Sum = SumMaker()

def istype(type_):
    return lambda arg: type(arg) == type_

def split(cond, list_):
    return filter(cond, list_), filter(~cond, list_)

class CoreClass(object):
    def __add__(self, *args):
        return Add(self, *args)
    
    def __pow__(self, arg):
        return Pow(self, arg)
    
    def eval(self, kind):
        obj = type(self).eval(self)
        obj = obj.eval()
        new_args = [ arg.eval(kind) for arg in obj.args]

    def __ior__(self, other):
        print( 'define', self, other)
                    
class ExprClass(CoreClass):
    def __init__(self, *args):
        self.args = args
        
    def __repr__(self):
        args = ', '.join( [ arg.__repr__() for arg in self.args] )
        return f"{self.__class__.__name__}({args})"
    
    def __hash__(self):
        return hash(repr(self))

    def __str__(self):
        return self.__repr__()
    
    def __eq__(self, other):
        print( "compare", self, other )
        if type(other) == AnyClass:
            return True
        if type(self) != type(other):
            return False
        
        return "not sure"
    
    def __rshift__(self, other):
        print( repr(self), ">>", repr(other) )
    
    def match(self, other):
        if isinstance(other, AnyClass):
            return [other, self]
        if type(self) == type(other):
            any_args, args = split(istype(AnyClass), other.args)
            for any_arg in any_args:
                any_, set(self.args) - set([ arg for arg in other.args if type(arg) != AnyClass])
            print( set(self.args) - set([ arg for arg in other.args if type(arg) != AnyClass]) )
        print( "no match", repr(other), ":", repr(self) )
        
#     def eval(self, kind=None):
        
        
class AnyClass(CoreClass):
    def __init__(self, name):
        self.name = name
        
    def __repr__(self):
        return f"Any({self.name})"
    
    def __eq__(self, other):
        return True
    
class SymbolClass(CoreClass):
    def __init__(self, name):
        self.name = name

    def __repr__(self):
        return f"Symbol({self.name})"

    def __str__(self):
        return fr"{self.name}"
        
class AddClass(ExprClass):
    def __str__(self):
        return ' + '.join( [ arg.__str__() for arg in self.args ] )
    
    def __eq__(self, other):
#         print( "Add compare", set(self.args), "and", set(other.args) )
        if set(self.args) - set(other.args) == set():
            return True
        return False
    
class SumClass(ExprClass):
    def __init__(self, lims, arg):
        super().__init__(lims, arg)
        self.lims = lims
        self.operand = arg
    
    def __str__(self):
        operand = fr"{ self.operand }"
        if type(self.operand) is AddClass:
            operand = delim(operand)
        return fr"\sum_{{ {self.lims} }}{{ {operand} }}"
    
class DeltaClass(ExprClass):
    pass