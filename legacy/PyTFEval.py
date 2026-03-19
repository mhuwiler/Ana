import numpy as np
import tensorflow as tf
import ctypes
import math
import h5py


class PyTFEval: 
    def __init__( self ):
        print('Initialising TFEvaluation')

        self.savedmodel = "../../data/batchsize_10/serialized"

        self.BATCHSIZE=10
        self.NUM_POINT = 20

        self.initCheck = False 

        self.debug = False

        #np.__config__.show()

        #os.environ['OPENBLAS_NUM_THREADS'] = '1'
        #os.environ['NUMEXPR_NUM_THREADS=1']

        self.Init()        

    def Initialise(self, model, batchsize, numpoints): 
        self.savedmodel = model
        self.BATCHSIZE = batchsize
        self.NUM_POINT = numpoints

        self.Init()

    def Init(self): 
        print("Using the following model: {}".format(self.savedmodel))
        self.sess = tf.Session(graph=tf.Graph())
        print("Inside session")

        tf.saved_model.loader.load(self.sess, [tf.saved_model.tag_constants.SERVING], self.savedmodel) #'../pretrained/{}'.format(FLAGS.name)
        self.output = self.sess.graph.get_tensor_by_name('Softmax_2:0') #Reshape_5:0
            
        #mock_data = np.ones((BATCHSIZE,NUM_POINT,NFEATURES),dtype=float)
        self.mock_label = np.ones((self.BATCHSIZE,self.NUM_POINT),dtype=float)
        #mock_glob = np.ones((self.BATCHSIZE,NGLOB),dtype=float)
        
    def pyArray (self, a):
        print ("Contents of a :")
        print (a)
        print(type(a))
        for b in a: 
            print(type(b))
        c = 0

        print(a.shape)

        return [1., 2., 3., 4., 5] #a.data_as(ctypes.POINTER(ctypes.c_double))

    def Eval(self, data): 
        #print("Evaluation called")

        #print(data.shape)

        data = np.transpose(data.reshape(13, 20))

        #print(data.shape)

        #print(data)

        data = np.expand_dims(data, 0)

        #print(data.shape)

        data[:, :, 2] = np.log(data[:, :, 2], out=np.zeros_like(data[:, :, 2]), where=(data[:, :, 2]!=0)) # taking log of pT 

        np.nan_to_num(data, copy=False)

        #self.CheckInput(data)

        return self.Evaluate(data)


    def Evaluate(self, data): 
        #print("In eval")
        shape = data.shape

        batchsize = shape[0]
        numpoints = shape[1]
        numvars = shape[2]

        #print(shape)

        assert(batchsize == 1), "ERROR: The fuction requires a single event!"
        if (numpoints > self.NUM_POINT): 
            data = data[:, 0:self.NUM_POINT, :]

        #print(shape)

        # update the shape after the consistency checks 
        shape = data.shape

        #print("{}, {}, {}".format(batchsize, numpoints, numvars))

        placeholder = np.zeros(shape)

        dim = self.BATCHSIZE - 1

        #print(placeholder)
        #print(placeholder.shape)

        placeholder = np.repeat(placeholder, dim, axis=0)

        #print(placeholder.shape)

        batch = np.concatenate((data, placeholder))

        #print(batch.shape)

        response = self.NN_response(batch)

        #with open("response.txt", "a") as outfile: 
        #    outfile.write(response)

        #response = [1., 2., 3., 4., 5]

        #print(response.shape)

        #print(response)

        result = response[0,:, :]

        #result = np.asarray(result, "d")

        if (self.debug): print(result)

        #print(result.shape)

        return result


    def NN_response(self, data):
        # the input to this definition is the first key in the h5 file
        # so far, this part is a blackbox for me
        
        # to distinguish if a track in short_data_set is found to be signal or background, read in the NN response (array of array of as many arrays as there are tracks. All of the event-arrays contain 3 probabilities, the first one tells you how likely it is a bg-track, the second if it's a signal-track and the third if it's a muon track.)  

        size = data.shape
        #print('input_size =', size)


        #print(data[0,1,:])
        
        
            
            
        
            
        pred_list = []
            
        for i in range(0, math.floor(len(data)/self.BATCHSIZE)): # TODO: remove
              
        # print("Evaluating from: {} to {}".format(i*self.BATCHSIZE+1, (i+1)*self.BATCHSIZE))
            feed_dict = {
                'Placeholder:0': data[i*self.BATCHSIZE:(i+1)*self.BATCHSIZE],
                'Placeholder_1:0': self.mock_label,
                'Placeholder_2:0': False,
            }
            '''
              if feed_dict['Placeholder:0'].shape == (10, 20, 13):
                pass
              else:
                print('feed_dict size of placeholder 0 different:', feed_dict['Placeholder:0'].shape)
            '''

        #print("Before running prediction")

        predictions = self.sess.run(self.output, feed_dict)

        #print("After running prediction")

        #print(predictions)
              
        #predictions store 1 value per particle. The shape is [1,100,2]. The probability as coming from the B decay is stored at position 2, for example, the probability for the first particle come from the B decay is stored at predictions[0,0,1], while for the second is stored at predictions[0,1,1] and so on.
        for pred in predictions:
            pred_list.append(pred)
            pred_array = np.asarray(pred_list, dtype=object)
        return np.asarray(pred_array, "d")

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


