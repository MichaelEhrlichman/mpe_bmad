import numpy as np
import ctypes as C
import scipy.optimize as soz

patch_maker = C.CDLL('/afs/psi.ch/project/slsbd/data/bmad/bmad_local/production/lib/libpatch_maker_4py.so')
funcprot = getattr(patch_maker,"patch_maker_4py")
funcprot.argtypes = [C.c_char_p, C.c_int, C.c_double, C.c_double, C.POINTER(C.c_double), C.POINTER(C.c_double)]

file = b"lat.bmad"
file_len = len(file)

def thefunc(xvec):
	field_scale = xvec[0]
	pitch = xvec[1]
	x = C.c_double(0.0)
	xp = C.c_double(0.0)
	funcprot(file, file_len, field_scale, pitch, x, xp)
	print(field_scale,pitch,x.value,xp.value)
	return x.value, xp.value

sol = soz.root(thefunc,[1.1, -0.037864398262])

print(sol)

