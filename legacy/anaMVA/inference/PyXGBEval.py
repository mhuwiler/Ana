from __future__ import division, print_function
import numpy as np
import pandas as pd
import sys
import ctypes
import math
import h5py
import xgboost as xgb
import pickle


class PyXGBEval: 
    def __init__( self ):
        print('Initialising XGBEvaluation')

        self.savedmodel = "../../data/batchsize_10/serialized"

        #classifier = pickle.load(open(options.model, "rb"))

        self.debug = False     

    def Initialise(self, model): 
        self.savedmodel = model

        self.Init()

    def Init(self): 
        print("Using the following model: {}".format(self.savedmodel))

        self.classifier = pickle.load(open(self.savedmodel, "rb"))

        sys.path.append('..')
        from Features import features_save

        self.features = [item[0] for item in features_save]
        

    def Eval(self, data): 
        #print("Evaluation called")

        dataforconvert = pd.Series(data)

        dataForEval = xgb.DMatrix(dataforconvert, feature_names=self.features)

        return self.Evaluate(dataForEval)


    def Evaluate(self, data): 

        #print(data)
        
        result = self.classifier.predict(data)

        #print(result)

        return result


    def CheckInput(self, dataframe): 
        # open the file and check if the dataframes are  identical
        print("Check")
        print(dataframe)
        if (self.initCheck != True): 
            self.InitCheck("../../data/firstAllTauFix.h5")
        
        print('data_set =', self.dataset.shape) 
        element = self.dataset[self.count]
        element = np.expand_dims(element, 0)
        print(element)
        self.count+=1
        if (not np.allclose(dataframe, element)): #assert(np.allclose(dataframe, self.dataset[self.count]))
            print("Warning: different values in arrays!")
        print(dataframe - element)
        assert(np.allclose(dataframe, element))
        print(dataframe.dtype)
        print(element.dtype)


    def InitCheck(self, filename): 
        self.file = h5py.File(filename, 'r') 
        self.dataset = self.file['data']
        self.initCheck = True
        self.count = 0


    """
    def EvaluateBatch(self, batch): 
        with tf.Session(graph=tf.Graph()) as sess:
            tf.saved_model.loader.load(sess, [tf.saved_model.tag_constants.SERVING], self.savedmodel) #'../pretrained/{}'.format(FLAGS.name)
            output = sess.graph.get_tensor_by_name('Softmax_2:0') #Reshape_5:0
            
            #mock_data = np.ones((BATCHSIZE,NUM_POINT,NFEATURES),dtype=float)
            mock_label = np.ones((self.BATCHSIZE,self.NUM_POINT),dtype=float)
            #mock_glob = np.ones((self.BATCHSIZE,NGLOB),dtype=float)

            feed_dict = {
                'Placeholder:0': batch,
                'Placeholder_1:0': mock_label,
                'Placeholder_2:0': False,
            }

            predictions = sess.run(output, feed_dict)

            return predictions
    """


