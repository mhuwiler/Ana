#!/usr/bin/env python

import PyTFEval
from PyTFEval import PyTFEval
import ROOT
import numpy


print("Testing TFEvaluation.py")


evaluation = PyTFEval()


# Test the function which takes a single row
import numpy as np
import h5py
print("Test placeholder")
f = h5py.File("../../data/firstAllTau.h5", 'r') #"/pnfs/psi.ch/cms/trivcat/store/user/mhuwiler/data/Analysis/taudnn/v10/test_TAU.h5"
maxnum = 100 #math.floor(len(f["data"])/10)*10
data_set = f['data'][0:maxnum]
print('data_set =', data_set.shape) 
label_set = f['label'][0:maxnum]
col = 1
#data_set[:,:,10] = data_set[:,:,11]

print(numpy.__version__)

canv1 = ROOT.TCanvas("canv1", "Canvas 1", 800, 600)
h1 = ROOT.TH1D("h1", "Event by event inference", 200, -0.1, 1.1)
# make predictions
for i in range(len(data_set[0:100])): 
    batch = data_set[i]
    label = label_set[i]
    #batch = np.expand_dims(batch, axis=0)
    print(batch.shape)

    pred_array = evaluation.Eval(batch)

    print(pred_array)

    column = pred_array[:,col]
    print(label.shape)
    column = column[label==1]
    dim = len(column)-1
    print("values {}".format(column))
    print(column.shape)

    #h1.FillN(dim, np.asarray(column, "d"), np.zeros(len(column)))
    for val in column: 
        h1.Fill(val)

h1.Draw("HIST")
canv1.Draw()
canv1.SetLogy()
canv1.Print("EventInference.pdf")


canv2 = ROOT.TCanvas("canv2", "Canvas 1", 800, 600)
h2 = ROOT.TH1D("h2", "Batch inference", 200, -0.1, 1.1)
# Test the legacy way 
for i in range(0, 90, 10): 
    batch = data_set[i:i+10]
    label = label_set[i:i+10]
    label = np.concatenate(label, axis=0)
    print(batch.shape)

    prediction = evaluation.NN_response(batch)

    print(prediction)

    print(prediction.shape)

    prediction = np.concatenate(prediction, axis=0)

    prediction = prediction[label==1]

    print(prediction.shape)

    column = prediction[:,col]
    dim = len(column)-1

    print("values {}".format(column))

    #h2.FillN(dim, np.asarray(column, "d"), np.zeros(len(column)))
    for var in column: 
        h2.Fill(var)

h2.Draw()
canv2.Draw("HIST")
canv2.SetLogy()
canv2.Print("BatchInference.pdf")


  