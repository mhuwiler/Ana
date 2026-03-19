#!/usr/bin/env python

# welcome to the file that is made of three individual python codes: read_pb_new.py, plot_eta.py and plot_pt.py
#from __future__ import division, print_function
import numpy as np
#import h5py
import tensorflow as tf
#from argparse import ArgumentParser
#import os, ast
#import sys
import math
#from sklearn.metrics import roc_curve
#import ROOT

print("Loading TFEvaluation.py")



class TFEvaluation: 
    def __init__( self ):
        print('Initialising TFEvaluation')

        self.savedmodel = "../../data/batchsize_10/serialized"

        self.BATCHSIZE=10
        self.NUM_POINT = 20

    def Condition(self, r):
        # Condition: returns a vector, where "True" matches the position of the rows with values 
        # False matches the positions of all-0-rows
        result = []
        for i in r: 
          z = not np.all(i == 0)
          result.append(z)
        return result

    def NN_response(self, data):
        # the input to this definition is the first key in the h5 file
        # so far, this part is a blackbox for me
        
        # to distinguish if a track in short_data_set is found to be signal or background, read in the NN response (array of array of as many arrays as there are tracks. All of the event-arrays contain 3 probabilities, the first one tells you how likely it is a bg-track, the second if it's a signal-track and the third if it's a muon track.)  

        size = data.shape
        print('input_size =', size)


        print(data[0,1,:])
        
        with tf.Session(graph=tf.Graph()) as sess:
            tf.saved_model.loader.load(sess, [tf.saved_model.tag_constants.SERVING], self.savedmodel) #'../pretrained/{}'.format(FLAGS.name)
            output = sess.graph.get_tensor_by_name('Softmax_2:0') #Reshape_5:0
            
            #mock_data = np.ones((BATCHSIZE,NUM_POINT,NFEATURES),dtype=float)
            mock_label = np.ones((self.BATCHSIZE,self.NUM_POINT),dtype=float)
            #mock_glob = np.ones((self.BATCHSIZE,NGLOB),dtype=float)
            
        
            
            pred_list = []
            
            latest_num = 0
            
            for i in range(0, math.floor(len(data)/self.BATCHSIZE)):
              
            # print("Evaluating from: {} to {}".format(i*self.BATCHSIZE+1, (i+1)*self.BATCHSIZE))
              feed_dict = {
                'Placeholder:0': data[i*self.BATCHSIZE:(i+1)*self.BATCHSIZE],
                'Placeholder_1:0': mock_label,
                'Placeholder_2:0': False,
              }
              
              latest_num = (i+1)*self.BATCHSIZE
              '''
              if feed_dict['Placeholder:0'].shape == (10, 20, 13):
                pass
              else:
                print('feed_dict size of placeholder 0 different:', feed_dict['Placeholder:0'].shape)
              '''
              predictions = sess.run(output, feed_dict)
              
            #predictions store 1 value per particle. The shape is [1,100,2]. The probability as coming from the B decay is stored at position 2, for example, the probability for the first particle come from the B decay is stored at predictions[0,0,1], while for the second is stored at predictions[0,1,1] and so on.
              for pred in predictions:
                pred_list.append(pred)
                pred_array = np.asarray(pred_list, dtype=object)
            return np.asarray(pred_array, "d")

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


    def Hello(self): 
        print("Hello!")
            
        

    def Evaluate(self, data): 
        shape = data.shape

        batchsize = shape[0]
        numpoints = shape[1]
        numvars = shape[2]

        assert(batchsize == 1), "ERROR: The fuction requires a single event!"
        if (numpoints > self.NUM_POINT): 
            data = data[:, 0:self.NUM_POINT, :]

        print(shape)

        # update the shape after the consistency checks 
        shape = data.shape

        print("{}, {}, {}".format(batchsize, numpoints, numvars))

        placeholder = np.zeros(shape)

        dim = self.BATCHSIZE - 1

        print(placeholder)
        print(placeholder.shape)

        placeholder = np.repeat(placeholder, dim, axis=0)

        print(placeholder.shape)

        batch = np.concatenate((data, placeholder))

        print(batch.shape)

        response = self.NN_response(batch)

        print(response.shape)

        return np.asarray(response[0,:], "d")


    def Eval(self, eta, phi, pt, q, DOCA2D, DOCA2DErr, DOCA, DOCAErr, dzToPV, dzToClosest, isAssociate, assocQuality): 
        df = np.zeros((1, 20, 13))

        df[1,:,1] = eta.Data()


        print(df)





    
    def plots_of_eta_and_pt(self, file_directory, modelpath):
        # for file_directory, the function requires a directory starting with / and ending with '.h5', the directory has to be written in quotation marks
        f = h5py.File(file_directory, 'r')
        maxnum = 100 #math.floor(len(f["data"])/10)*10
        data_set = f['data'][0:maxnum]
        print('data_set =', data_set.shape) 
        label_set = f['label'][0:maxnum]
        
        # here is the important difference: the concatenated data_set and label_set
        data_set_cat = np.concatenate(data_set, axis=0)
        label_set_cat = np.concatenate(label_set, axis=0)
        
        # now, we don't need a loop anymore to take out the values that are interesting for us => all rows that are not completely 0
        cond_vector = self.Condition(data_set_cat)
        short_data_set = data_set_cat[cond_vector]
        short_label_set = label_set_cat[cond_vector]
        print('short label set =', short_label_set.shape)
        
        # read in the pt values: short_data_set[ : ,[2]], until now, we have log(pt), but we now want to have pt
        short_data_set[ : ,[2]] = np.exp(short_data_set[ : ,[2]])
        
        # make predictions
        pred_array = self.NN_response(data_set)
        
        # concatenate again, this time the pred_array, to put into shape (batchsize*20):
        pred_array_cat = np.concatenate(pred_array, axis=0)
        
        # also here: take only the values that are interesting
        short_pred_array = np.squeeze(pred_array_cat[cond_vector])
        print('short pred array =', short_pred_array.shape)
        
        #return(background, signal, pred_background, pred_signal, short_label_set, short_pred_array)
        return [short_label_set, short_pred_array, short_data_set]


    def Test(self): 
        #background_test, signal_test, bkg_pred_test, signal_pred_test, labels_test, predictions_test = plots_of_eta_and_pt('/work/alorenze/bachelorthesis/TAU/test_TAU.h5', 
        #"/work/alorenze/bachelorthesis/taudnn/logs/fourthtraining/serialized/")
        resp = self.plots_of_eta_and_pt('../../data/v10/test_TAU.h5', 
        "../../data/batchsize_10/serialized")
        label_test = resp[0]
        predictions_test = resp[1]
        dataset = resp[2]
        print('shape of labels =', label_test.shape)
        print('shape of preds =', predictions_test.shape)
        print('predictions =', predictions_test) 
        signal_mask = label_test==1
        predictions_test = predictions_test[signal_mask]
        predictions_test = predictions_test.tolist()
        label_test = label_test.tolist()

        #print(predictions_test.shape)

        # ROOT plot 
        """
        canvas = ROOT.TCanvas("canvas", "canvas", 800, 600)
        histo = ROOT.TH1D("histo", "histo", 200, -0.1, 1.1)

        #histo.FillN(len(predictions_test)-1, np.asarray(predictions_test, "d"), np.zeros(len(predictions_test)))
        for values in predictions_test: 
            print(values)
            #histo.Fill(val)

            #histo.FillN(len(val)-1, np.asarray(val, "d"), np.zeros(len(val)))
            for val in values: 
                histo.Fill(val)

        histo.Draw("HIST")
        canvas.Draw()
        canvas.Print("SignalDistribution.pdf")

        with h5py.File('TestResponse.h5', "w") as f: 
            f.create_dataset('predictions', data=predictions_test)
            f.create_dataset('label', data=label_test)
        """



  