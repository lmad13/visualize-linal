import numpy as np
from mpl_toolkits.mplot3d import axes3d
import matplotlib.pyplot as plt
from System3D import System
import sys
print 'Beginning Testing of Visulalization'

if len(sys.argv)<2:
    print 'Improper usage! \nProper usage: python VisualizeLinearAlgebra.py testInput.data'
    print 'Exiting!'
else:
    fileName=sys.argv[1]
    #For the tetramer
    bonding=[(0,4),(0,5),(0,12),(1,6),(1,7),(1,10),(2,8),(2,9),(2,11),(3,12),(3,10),(3,11)]
    system=System()
    NameArray,posArray=system.loadMolecules(fileName)
    for Names, pos in zip(NameArray, posArray):
        system.addMolecule(Names,pos,bonding)
        system.addHCartDispAxesVects(pos,'tetramer')
        print 'added H cart Displacement'
        system.display()
        system.resetSys()
    
