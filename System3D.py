#This class is helpful for visualizing the linear algebra operations that I'm performing on my walkers. This might be usable for both the trimer and tetramer.


import numpy as np
from mpl_toolkits.mplot3d import axes3d
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
import glob


class System(object):
#
# The System object can be thought of as 3D cartesian space with defined points and vectors.
# The System can be added to or visualized.  After Visualization, the matplotlib plt is automatically reset. 
#
    def __init__(self):
        # A System is a set of atoms, bonds, and vectors that can be rendered with matplotlib.pyplot's axes3d library
        self.fig = plt.figure()
        self.ax=self.fig.gca(projection='3d')
        
        #So that the atoms appear in reasonable colors
        self.colordict={}
        self.colordict['C']='k' #carbon is black
        self.colordict['H']='w' #hydrogen white
        self.colordict['O']='r' #oxygen red, of course

        self.functAddList=[] #only rendering commands can go in here
        self.paramAddList=[] #parameters for the rendering commands go here

    
    #Functions for adding of elements to the system
    def addAtom(self,Name,x,y,z):
        self.functAddList.append(self.renderAtom)
        self.paramAddList.append([Name,x,y,z])        

    def addMolecule(self,Names,pos,bonding):
        for n,coord in zip(Names,pos):
            self.addAtom(n,coord[0],coord[1],coord[2])
        print 'added Molecule'
        verts=[]
        zs=[]
        for pair in bonding:
            self.addBond(pos[pair[0]],pos[pair[1]])

    def addBond(self,atomStart, atomEnd):
        self.functAddList.append(self.renderBond)
        self.paramAddList.append([atomStart, atomEnd])

    def addHCartDispAxesVects(self,pos,moleculeName):
        self.functAddList.append(self.renderCartDispAxesVects)
        self.paramAddList.append([pos,moleculeName])


    #Functions for rendering the atoms/bonds/vectors
    def renderAtom(self,args):
        [Name, x, y, z] = args
        self.ax.scatter(x,y,z,c=self.colordict[Name], s=100)
    def renderBond(self,args):
        [atomStart,atomEnd]=args
        self.ax.plot([atomStart[0],atomEnd[0]],[atomStart[1],atomEnd[1]],[atomStart[2],atomEnd[2]],color='k')

    
    def renderCartDispAxesVects(self,args):
        #adds xyz coordinates at the center point (c) between solvating a O and the central O 
        #adds a black displacement vector pointing between c and the shared proton
        #adds thicker xyz component vectors corresponding to the black displacement vector

        [pos,moleculeName]=args
        if moleculeName=='tetramer':
            #If the order of your atoms is different than the Madison Version, the reordering should happen here.
            O1=pos[0] 
            O2=pos[1]
            O3=pos[2]
            OC=pos[3]  #central oxygen 
            H1=pos[12]  #omg
            H2=pos[10]
            H3=pos[11]
                    
            c1=(O1+OC)/2.0  #center points
            c2=(O2+OC)/2.0
            c3=(O3+OC)/2.0
            
            X1prime=O1-OC  #vectors along OO bonds
            X1prime=X1prime/np.linalg.norm(X1prime)
            X2prime=O2-OC
            X2prime=X2prime/np.linalg.norm(X2prime)
            X3prime=O3-OC
            X3prime=X3prime/np.linalg.norm(X3prime)
            
            normal1=np.cross(O1-OC,O3-OC) #normal to the O_c-O_a-O_b planes
            normal2=np.cross(O2-OC,O1-OC)
            normal3=np.cross(O3-OC,O2-OC)
            
            Z1prime=normal1/np.linalg.norm(normal1)
            Z2prime=normal2/np.linalg.norm(normal2)
            Z3prime=normal3/np.linalg.norm(normal3)
            Y1prime=np.cross(Z1prime,X1prime)  #in the plane of the O_c-O_b-O_a planes
            Y2prime=np.cross(Z2prime,X2prime)
            Y3prime=np.cross(Z3prime,X3prime)

            Y1prime=Y1prime/np.linalg.norm(Y1prime)
            Y2prime=Y2prime/np.linalg.norm(Y2prime)
            Y3prime=Y3prime/np.linalg.norm(Y3prime)
            
            axes=[[X1prime,Y1prime,Z1prime],[X2prime,Y2prime,Z2prime],[X3prime,Y3prime,Z3prime]]
            #rending reference axes
            for atomC,axis in zip([c1,c2,c3],axes):
                for direction,c in zip(axis,['r','g','b']):
                    start=atomC
                    end=start+direction
                    coords=zip(start,end)
                    self.ax.plot(coords[0],coords[1],coords[2], linewidth=1.0,color=c)
            
            #calculating shared proton displacements
            H1dispx=np.dot(X1prime,H1-c1)
            H1dispy=np.dot(Y1prime,H1-c1)
            H1dispz=np.dot(Z1prime,H1-c1)
            
            H2dispx=np.dot(X2prime,H2-c2)
            H2dispy=np.dot(Y2prime,H2-c2)
            H2dispz=np.dot(Z2prime,H2-c2)
            
            H3dispx=np.dot(X3prime,H3-c3)
            H3dispy=np.dot(Y3prime,H3-c3)
            H3dispz=np.dot(Z3prime,H3-c3)

            #rendering displacement
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


    #Functions related to displaying or resetting the plotting environment

    def display(self):
        #executes the plotting commands and displays
        for f,p in zip(self.functAddList,self.paramAddList):
            f(p)
        plt.show()
        print 'Displayed the System!'


    def resetFig(self):
        #Resets the plotting environment
        self.fig = plt.figure()
        self.ax=self.fig.gca(projection='3d')

    def resetSys(self):
        #Resets the plotting environment AND the list of atoms/bond/vectors 
        self.fig = plt.figure()
        self.ax=self.fig.gca(projection='3d')
        self.functAddList=[]
        self.paramAddList=[]
        print 'reset Figure'

    #Funtions for loading in from xyz files
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
        names=np.array(names).reshape((10,nAtoms))
        return names, coords
