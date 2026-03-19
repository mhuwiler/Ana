import numpy as np
import tensorflow as tf
import ctypes


class TestPyTFEval: 
    def __init__( self ):
        print('Initialising TFEvaluation')

        self.savedmodel = "../../data/batchsize_10/serialized"

        self.BATCHSIZE=10
        self.NUM_POINT = 20
        
    def pyArray (self, a):
    	print ("Contents of a :")
    	print (a)
    	print(type(a))
    	for b in a: 
    		print(type(b))
    	c = 0
    	return a #a.data_as(ctypes.POINTER(ctypes.c_double))
