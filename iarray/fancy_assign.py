Q = {}

class Symbol:
    global Q
    def __init__(self, key ):
        self.key = key
        
    def __getitem__(self, ind):
        key = f"{self.key}_{{ {ind} }}"
        return Symbol(key)
    
    def __lshift__(self, value):
        Q[str(self.key)] = value
        
    def __repr__(self):
        return self.key
                
def diff( q, *args ):
    num = fr"{{\rm d}}\! {q.key}"
    den = [ fr"{{\rm d}}\! {arg.key}" for arg in args ]
    key = fr"\frac{{ {num} }}{{ {''.join(den)} }}"
    return Symbol(key)