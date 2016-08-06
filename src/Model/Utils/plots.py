import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np
from numpy import genfromtxt
import matplotlib.cbook as cbook

data = genfromtxt('../../Data/twod.csv', delimiter=',')
classes = genfromtxt('../../Data/simulation.csv')
classesP = genfromtxt('../../Data/simulationP.csv')
print(classesP)

plt.scatter(data[:,0], data[:,1], c=classes, marker='o')
plt.show(block='False')
plt.figure()
plt.scatter(data[:,0], data[:,1], c=classesP, marker='o')
plt.show()