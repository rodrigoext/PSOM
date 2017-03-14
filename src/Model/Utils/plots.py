import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.cm as cm
import numpy as np
from numpy import genfromtxt
import matplotlib.cbook as cbook

has_data = False
#data = genfromtxt('../../Data/data.csv', delimiter=',')
neurons = genfromtxt('../../Data/codebook.csv', delimiter=',')
ustar = genfromtxt('../../Data/ustar.csv', delimiter=',')
um = genfromtxt('../../Data/um.csv', delimiter=',')
ustarw = genfromtxt('../../Data/ustarw.csv', delimiter=',')
immersion = genfromtxt('../../Data/immersion.csv', delimiter=',')
pmatrix = genfromtxt('../../Data/pmatrix.csv', delimiter=',')
#classes = genfromtxt('../../Data/simulation.csv')
#classesP = genfromtxt('../../Data/simulationP.csv')
classesN = genfromtxt('../../Data/classes.csv')
classesI = genfromtxt('../../Data/classesIm.csv')
#print(classesP)

if has_data :
	data = genfromtxt('../../Data/data.csv', delimiter=',')
	classesP = genfromtxt('../../Data/simulationP.csv')
	plt.scatter(data[:,0], data[:,1], c=classesP, marker='o')
	plt.draw()
plt.figure()
plt.scatter(neurons[:,0], neurons[:,1], c=classesN, marker='o')
plt.draw()
#fig = plt.figure()
#ax = Axes3D(fig)
#ax.scatter(neurons[:,0], neurons[:,1], neurons[:,2], 'z', 60, c=classesN, #marker='o')
plt.draw()
plt.figure()
plt.scatter(neurons[:,0], neurons[:,1], c=classesI, marker='o')
plt.draw()
plt.figure()
plt.imshow(um)
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
