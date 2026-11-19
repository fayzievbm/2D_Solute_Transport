"""Process/thread-local IEEE underflow mode; restored on exit.

Only subnormal double values (<2.225074e-308 in magnitude) are flushed.
No concentration threshold, negativity clipping, or fastmath is used.
Microsoft documentation: https://learn.microsoft.com/en-us/cpp/c-runtime-library/reference/controlfp-s
"""
from contextlib import contextmanager
import ctypes,sys
import numpy as np
NORMAL_MIN=float(np.finfo(np.float64).tiny)

@contextmanager
def flush_subnormals():
    if sys.platform!='win32':
        yield {'enabled':False,'normal_min':NORMAL_MIN};return
    lib=ctypes.CDLL('ucrtbase');fn=lib._controlfp_s
    fn.argtypes=[ctypes.POINTER(ctypes.c_uint),ctypes.c_uint,ctypes.c_uint];fn.restype=ctypes.c_int
    old=ctypes.c_uint();assert fn(ctypes.byref(old),0,0)==0
    new=ctypes.c_uint();assert fn(ctypes.byref(new),0x01000000,0x03000000)==0
    try:yield {'enabled':True,'normal_min':NORMAL_MIN,'control_word':new.value,'scope':'current simulation thread only'}
    finally:
        result=ctypes.c_uint();assert fn(ctypes.byref(result),old.value,0x03000000)==0
