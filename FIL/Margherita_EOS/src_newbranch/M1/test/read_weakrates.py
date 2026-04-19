import numpy as np
from enum import Enum

class Header(Enum):
    TEMP=0
    Q = 1
    K_S = 2
    K_A = 3
    EPS2_FLUID = 4
    EPS2_LEAK  = 5
    EPS  = 6
    ETA  = 7
    ETAP = 8
    ETAN = 9
    FTAU = 10
    TAU  = 11
    
def read_weakrates(fname):
    raw_data = np.loadtxt(fname)
    out = dict()
    for name,member in Header.__members__.items():
        out[name] = raw_data[:,member.value]
    return out 
