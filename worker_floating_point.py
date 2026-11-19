"""Apply the verified subnormal mode to each simulation worker, then restore it."""
import ctypes,sys
from contextlib import contextmanager
import numpy as np
from numba import njit,prange,get_num_threads,get_thread_id
if sys.platform=='win32':
    _lib=ctypes.CDLL('ucrtbase');_control=_lib._controlfp
    _control.argtypes=[ctypes.c_uint,ctypes.c_uint];_control.restype=ctypes.c_uint

    @njit(parallel=True)
    def _set(saved):
        for k in prange(4*len(saved)):
            tid=get_thread_id()
            if saved[tid]==0xffffffff:saved[tid]=_control(0,0)
            _control(0x01000000,0x03000000)

    @njit(parallel=True)
    def _restore(saved):
        for k in prange(4*len(saved)):
            tid=get_thread_id();_control(saved[tid]&0x03000000,0x03000000)

@contextmanager
def worker_subnormals():
    if sys.platform!='win32':yield {'enabled':False};return
    saved=np.full(get_num_threads(),0xffffffff,dtype=np.uint32)
    _set(saved)
    assert (saved!=0xffffffff).all(),'Not every worker was configured'
    try:yield {'enabled':True,'workers':len(saved),'saved_control_words':[int(x) for x in saved],'normal_min':2.2250738585072014e-308,'scope':'simulation worker threads; restored on exit'}
    finally:_restore(saved)
