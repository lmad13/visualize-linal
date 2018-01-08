import numpy as np
from mpl_toolkits.mplot3d import axes3d
import matplotlib.pyplot as plt
from System3D import System
import sys
print 'Beginning Visulalization'



pos=[ [ -2.6904576112310954      ,  3.9500841728884906      , 0.75942179065411275      ],
 [  -1.7952128568450187     ,  -4.2353599559094803     ,   2.2411643264099834E-002],
 [   4.4422336549277324     ,  0.45406749171931748     ,  -3.8829021991088897E-002],
 [   4.8201545007033886E-002,  -8.7170580484709465E-002, -0.70059224887304339     ],
 [  -2.8171347706359655     ,   4.9583419999527107     , -0.63003976536089690     ],
 [  -3.4922677270157867     ,   4.7095394465380611     ,   2.4931723882396208     ],
 [  -1.7642989702093939     ,  -5.4646609736834870     ,  -1.5007325569448491     ],
 [  -2.1195419108788669     ,  -5.0911480361002646     ,   1.4918834096384197     ],
 [   6.0223081596814829     , -0.13357863218486773     ,  0.44438575332504177     ],
 [   4.4444214987729662     ,  0.12863016588281256     ,  -1.7316245025870030     ],
 [  -1.0072962960941543     ,  -1.7778651133095729     , -0.12200658346766519     ],
 [   1.9097746075432118     ,  0.13049311309059597     ,  -3.3865528683480473E-002],
 [  -1.2515843476459452     ,   1.2448614220389171     ,  -1.0842844570050636     ]]

Names='O O O O H H H H H H H H H'
bonding=[(0,4),(0,5),(0,12),(1,6),(1,7),(1,10),(2,8),(2,9),(2,11),(3,12),(3,10),(3,11)]
pos=np.array(pos)
Names=Names.split()
fileName=sys.argv[1]
system=System()
NameArray,posArray=system.loadMolecules(fileName)
for Names, pos in zip(NameArray, posArray):
    system.addMolecule(Names,pos,bonding)
    system.addHCartDispAxesVects(pos,'tetramer')
    print 'added H cart Displacement'
    system.display()
    system.resetFig()


