from matplotlib import pyplot as plt
import numpy as np
import argparse

epochs = np.load('output.epochs.npy')
freqs = np.load('output.freqs.npy')
logpost = np.load('output.post.npy')


import csv

with open("logposts.csv","w+") as my_csv:
    csvWriter = csv.writer(my_csv,delimiter=',')
    csvWriter.writerows(logpost)

epochs.tofile('epochs.csv',sep=',',format='%10.5f')