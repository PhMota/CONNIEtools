from .rule import Rule, RuleList

def translate( arg ):
    if isinstance(arg, Operable):
        return arg
    if isinstance(arg, (int, float)):
        if arg < 0:
            return Neg(Number(abs(arg)))
        return Number(arg)
    if isinstance(arg, str):
        return Symbol(arg)
    raise Exception('translate not implemented for', arg.__class__)

class Add:
    def __init__(self, label):
        
    def __add__(self, arg):
        return Add(self, arg)

_instances = {}
class Operable:
    '''
    this is the basic class that defines the interface with Python
    mainly, it converts Python operators into Functions and inputs
    '''
    def __hash__(self):
        return hash(self.label)
    
    def __new__(cls, label):
        _hash = hash(label)
        if not _hash in _instances.keys():
            new_obj = object().__new__(cls)
            new_obj.label = label
            new_obj.rules = RuleList()
            _instances[_hash] = new_obj
            return new_obj
        return _instances[_hash]

    def __init__(self, label):
        pass
    
    @classmethod
    def instance(cls, obj):
        return isinstance(obj, cls)
    
    @classmethod
    def list(cls):
        print( cls._instances )
    
    def __del__(self):
        print('destroyed obj', str(self), hash(self.label) )
        del _instances[hash(self.label)]
        del self
        
    def __repr__(self):
        return str(self.label)

    def __str__(self):
        return str(self.name)
        
    def __eq__(self, arg):
        return hash(self.label) == hash(translate(arg).label)
    
    def __neg__(self):
        from .neg import Neg
        return Neg( self )
    
    def __add__(self, arg):
        return Add( self, translate(arg) )

    def __radd__(self, arg):
        return Add( translate(arg), self)

    def __sub__(self, arg):
        return Add( self, Neg(translate(arg)) ) 

    def __rsub__(self, arg):
        return Add( translate(arg), Neg(self) )
    
    def __mul__(self, arg):
        return Mul( self, translate(arg) )

    def __rmul__(self, arg):
        return Mul( translate(arg), self )

    def __truediv__(self, arg):
        return Mul( self, Pow(translate(arg), translate(-1)) )

    def __rtruediv__(self, arg):
        return Mul( translate(arg), Pow(self, translate(-1)) )

    def __pow__(self, arg):
        return Pow( self, translate(arg))
    
    def __rpow__(self, arg):
        return Pow( translate(arg), self )
    
    def __call__(self, *args):
        if hasattr(args[0], '__iter__'):
            args = args[0]
        return AppliedFunction( self, *map(translate, args) )
    
    def __getitem__(self, indices):
        return Indexed(self, indices)

    def eval(self, tab=0, verbose=False):
        obj = self
        while True:
            new_obj = Operable.recursive_eval(obj, tab=0, verbose=verbose)
            if new_obj == obj:
                break
            obj = new_obj
        return new_obj
    
    @classmethod
    def recursive_eval(cls, obj, tab=0, verbose=False):
        new_obj = obj.rules( obj, verbose=verbose )
        new_obj = new_obj.__class__.rules( new_obj, verbose=verbose )
        if hasattr(new_obj, 'args'):
            new_args = []
            for arg in new_obj.args:
                arg = Operable.recursive_eval(arg, tab=tab+1, verbose=verbose)
                new_args.append( arg )
            new_obj = new_obj.operator(*new_args)
        return new_obj

class Expr(Operable):
    def __new__(cls, *args):
        label = fr'{cls.__name__}({", ".join(map(repr, args))})'
        obj = super().__new__(cls, label)
        if not hasattr(obj, 'operator'):
            obj.operator = cls
            obj.args = args
            obj.rules = RuleList()
        return obj
    
    def __init__(self, *args, **kwargs):
        pass
    
    def __str__(self):
        return self.label


def rules_mul(obj, verbose=False):
    args = flatten_args(obj)    
    terms = Counter()
    numeric = 1
    for arg in args:
        if isinstance(arg, Add):
            return Add(*[ Mul(add_arg, *[ other for other in args if not other == arg ]) for add_arg in arg.args ])
        if isinstance(arg, Neg):
            iarg = arg.args[0]
            numeric *= -1
            if isinstance(iarg, Number):
                numeric *= iarg.value
            else:
                terms.update({iarg:1})
        elif isinstance(arg, Number):
            numeric *= arg.value
        elif isinstance(arg, Pow):
            if isinstance(arg.power, Neg):
                if isinstance(arg.power.arg, Number):
                    terms.update({arg.base: -arg.power.arg.value})
            elif isinstance(arg.power, Number):
                terms.update({arg.base: arg.power.value})
            else:
                terms.update({arg.base: arg.power})
        else:
            terms.update({arg:1})
    new_args = [ s**v if v != 1 else s for s, v in terms.items() if v != 0]
    if numeric == 0:
        return Zero
    if len(new_args) == 0:
        return Number(numeric)
    if len(args) == 1 and numeric == 1:
        return new_args[0]
    if len(args) == 1 and numeric == -1:
        return Neg(new_args[0])
    if numeric == 1:
        return Mul( *new_args )
    if numeric == -1:
        return Neg(Mul( *new_args ))
    ret = Mul( Number(numeric), *new_args )
    return ret

class Mul(Expr):
    @staticmethod
    def add_order():
        return Symbol.add_order()+2
    @staticmethod
    def mul_order():
        return Symbol.mul_order()+1
    
    rules = RuleList(
        Rule(
            r'mul rules',
            lambda obj: True,
            lambda obj, v: rules_mul(obj, v)
        ),
    )

    def __new__(cls, *args, **kwargs):
        if len(args) == 0:
            return Zero
        if len(args) == 1:
            return args[0]
        if Zero in args:
            return Zero
        if One in args:
            return Mul(*[arg for arg in args if not arg == One])
        args = sorted( args, key=lambda x: x.mul_order() )
        obj = super().__new__(cls, *args )
        return obj
    
    def __init__(self, *args, **kwargs):
        pass
    
    def __str__(self):
        s = ''
        for arg in self.args:
            if isinstance(arg, (Add, Neg, Mul)):
                s += delim(str(arg))
            else:
                s += str(arg)
        return s
