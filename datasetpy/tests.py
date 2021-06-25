import numpy as np
from uncertainties import unumpy as un

from .core import *

def test_is_str():
    assert is_str(" ")
    assert not is_str([" "])
    
def test_is_list_or_tuple():
    assert is_list_or_tuple(["a"])
    assert is_list_or_tuple(("a",))
    assert not is_list_or_tuple("a")
    assert not is_list_or_tuple({"a": "b"})
    
def test_has_std_devs():
    assert not has_std_devs( np.arange(0,1,.1) )
    assert not has_std_devs( un.uarray(np.arange(0,1,.1), 0) )
    assert has_std_devs( un.uarray(np.arange(0,1,.1), .01) )

def  test_Coord__init__():
    c_np = Coord( 
        name = "my_coord", 
        data = np.arange(0,1,.1),
        dims = "my_dim"
    )
    c_np0 = Coord( 
        name = "my_coord", 
        data = np.arange(0,1,.1),
        err = 0,
        dims = "my_dim"
    )
    c_np1 = Coord( 
        name = "my_coord", 
        data = np.arange(0,1,.1),
        err = 1,
        dims = "my_dim"
    )
    c_un0 = Coord( 
        name = "my_coord", 
        data = un.uarray( np.arange(0,1,.1), 0 ),
        dims = "my_dim"
    )
    c_un1 = Coord( 
        name = "my_coord", 
        data = un.uarray( np.arange(0,1,.1), 1 ),
        dims = "my_dim"
    )
    assert str(c_np) == str(c_un0)
    assert str(c_np1) == str(c_un1)
    assert np.all( c_np == c_un0 )