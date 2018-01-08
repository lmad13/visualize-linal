#This class is helpful for visualizing the linear algebra operations that I'm performing on my walkers. This might be usable for both the trimer and tetramer.


import numpy as np
from mpl_toolkits.mplot3d import axes3d
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
import glob
class System(object):
    def __init__(self):
        self.fig = plt.figure()
        self.ax=self.fig.gca(projection='3d')
        self.colordict={}
        self.colordict['C']='k'
        self.colordict['H']='w'
        self.colordict['O']='r'
    def resetFig(self):
        self.fig = plt.figure()
        self.ax=self.fig.gca(projection='3d')
        print 'reset Figure'
    def addAtom(self,Name,x,y,z):
        self.ax.scatter(x,y,z,c=self.colordict[Name], s=100)
    def addBond(self,atomStart, atomEnd):
        self.ax.plot([atomStart[0],atomEnd[0]],[atomStart[1],atomEnd[1]],[atomStart[2],atomEnd[2]],color='k')
    def display(self):
        plt.show()
        print 'displayed!'
    def addMolecule(self,Names,pos,bonding):
        for n,coord in zip(Names,pos):
            self.addAtom(n,coord[0],coord[1],coord[2])
        print 'added Molecule'
        verts=[]
        zs=[]
        for pair in bonding:
            self.addBond(pos[pair[0]],pos[pair[1]])
    def loadMolecules(self,fileName):
        #grabs 10 Molecules
        f=open(fileName, 'r')
        line= f.readline()
        nAtoms=int(line.split()[0])
        print 'there are ',nAtoms
        rptunt=nAtoms+4
        coords=np.zeros((10,nAtoms,3))
        names=[]

        f=open(fileName, 'r')
        for j in range(10):
            f.readline()# Natoms
            f.readline() #Comments
            for a in range(nAtoms):
                data=f.readline().split()
                names.append(data[0])

                coords[j,a,:]=[float(data[1]),float(data[2]), float(data[3])]
            f.readline() #blank line
        print coords.shape
        names=np.array(names).reshape((10,nAtoms))
        print names.shape
        return names, coords

    def addHCartDispAxesVects(self,pos,moleculeName):
        if moleculeName=='tetramer':
            O1=pos[0]
            O2=pos[1]
            O3=pos[2]
            OC=pos[3]  #central oxygen 
            H1=pos[12]  #omg
            H2=pos[10]
            H3=pos[11]
                    
            c1=(O1+OC)/2.0
            c2=(O2+OC)/2.0
            c3=(O3+OC)/2.0
            
            X1prime=O1-OC
            X1prime=X1prime/np.linalg.norm(X1prime)
            X2prime=O2-OC
            X2prime=X2prime/np.linalg.norm(X2prime)
            X3prime=O3-OC
            X3prime=X3prime/np.linalg.norm(X3prime)
            
            #normal1=np.cross(O1-OC,O2-OC)
            #normal2=np.cross(O2-OC,O3-OC)
            #normal3=np.cross(O3-OC,O1-OC)

            normal1=np.cross(O1-OC,O3-OC)
            normal2=np.cross(O2-OC,O1-OC)
            normal3=np.cross(O3-OC,O2-OC)
            
            Z1prime=normal1/np.linalg.norm(normal1)
            Z2prime=normal2/np.linalg.norm(normal2)
            Z3prime=normal3/np.linalg.norm(normal3)
            Y1prime=np.cross(Z1prime,X1prime)
            Y2prime=np.cross(Z2prime,X2prime)
            Y3prime=np.cross(Z3prime,X3prime)

            Y1prime=Y1prime/np.linalg.norm(Y1prime)
            Y2prime=Y2prime/np.linalg.norm(Y2prime)
            Y3prime=Y3prime/np.linalg.norm(Y3prime)
            


            axes=[[X1prime,Y1prime,Z1prime],[X2prime,Y2prime,Z2prime],[X3prime,Y3prime,Z3prime]]
            for atomC,axis in zip([c1,c2,c3],axes):
                for direction,c in zip(axis,['r','g','b']):
                    start=atomC
                    end=start+direction
                    coords=zip(start,end)
                    self.ax.plot(coords[0],coords[1],coords[2], linewidth=1.0,color=c)
            

                     
                     
            H1dispx=np.dot(X1prime,H1-c1)
            H1dispy=np.dot(Y1prime,H1-c1)
            H1dispz=np.dot(Z1prime,H1-c1)
            
            H2dispx=np.dot(X2prime,H2-c2)
            H2dispy=np.dot(Y2prime,H2-c2)
            H2dispz=np.dot(Z2prime,H2-c2)
            
            H3dispx=np.dot(X3prime,H3-c3)
            H3dispy=np.dot(Y3prime,H3-c3)
            H3dispz=np.dot(Z3prime,H3-c3)

            axes=[[H1dispx*X1prime,H1dispy*Y1prime,H1dispz*Z1prime],
                  [H2dispx*X2prime,H2dispy*Y2prime,H2dispz*Z2prime],
                  [H3dispx*X3prime,H3dispy*Y3prime,H3dispz*Z3prime]]
            for atomC,axis in zip([c1,c2,c3],axes):
                for direction,c in zip(axis,['r','g','b']):
                    start=atomC
                    end=start+direction
                    coords=zip(start,end)
                    self.ax.plot(coords[0],coords[1],coords[2], linewidth=3.0,color=c)
            for atomC,atomH in zip([c1,c2,c3],[H1,H2,H3]):
                coords=zip(atomC,atomH)
                self.ax.plot(coords[0],coords[1],coords[2], linewidth=3.0,color='k')
