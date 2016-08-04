import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np
from numpy import genfromtxt
import matplotlib.cbook as cbook

data = genfromtxt('../../Data/twod.csv', delimiter=',')
classes = genfromtxt('../../Data/simulation.csv')
print(classes)

plt.scatter(data[:,0], data[:,1], c=classes, marker='o')

plt.show()