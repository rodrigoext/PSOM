import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np
from numpy import genfromtxt
import matplotlib.cbook as cbook

data = genfromtxt('../../Data/twod.csv', delimiter=',')
neurons = genfromtxt('../../Data/codebook.csv', delimiter=',')
ustar = genfromtxt('../../Data/ustar.csv', delimiter=',')
ustarw = genfromtxt('../../Data/ustarw.csv', delimiter=',')
immersion = genfromtxt('../../Data/immersion.csv', delimiter=',')
pmatrix = genfromtxt('../../Data/pmatrix.csv', delimiter=',')
#classes = genfromtxt('../../Data/simulation.csv')
classesP = genfromtxt('../../Data/simulationP.csv')
print(classesP)

# plt.scatter(data[:,0], data[:,1], c=classes, marker='o')
# plt.show(block='False')
# plt.figure()
plt.scatter(data[:,0], data[:,1], c=classesP, marker='o')
plt.draw()
plt.figure()
plt.scatter(neurons[:,0], neurons[:,1], marker='*')
plt.draw()
plt.figure()
plt.imshow(pmatrix)
plt.draw()
plt.figure()
plt.imshow(ustar)
plt.draw()
plt.figure()
plt.imshow(immersion)
plt.draw()
plt.figure()
plt.imshow(ustarw)
plt.show()
