

import numpy as np
from numpy.random import *

seed()


class Ising():
    def __init__(self, L, temp):
        self.L = L
        self.T = temp

    def Initialize(self):
        self.state = np.zeros((self.L, self.L))
        for i in range(self.L):
            for j in range(self.L):
                if randint(0, 2) > 0.5:
                    self.state[i][j] = 1  # Dipole has spin up
                else:
                    self.state[i][j] = -1
        return self.state
        
    def Hamiltonian(self):
        self.E = np.zeros((self.L, self.L))
        for i in range(self.L):
            for j in range(self.L):
                self.E[i, j] = -self.state[i, j]*(self.state[(i+1)%self.L, j] + self.state[i,(j+1)%self.L] + self.state[(i-1)%self.L, j] + self.state[i,(j-1)%self.L])

        #Taking care of spins except edges and corners i.e.  center = center*(top+bottom+left+right)
        #self.E[1:-1,1:-1] = -self.state[1:-1,1:-1]*(self.state[2:,1:-1]+self.state[:-2,1:-1]+self.state[1:-1,:-2]+self.state[1:-1,2:])
		#Taking care of edges
		# #Left
		# self.E[1:-1,0] = -self.state[1:-1,0]*(self.state[2:,-1]+self.state[:-2,-1]+self.state[1:-1,-1]+self.state[1:-1,2])
		# #right
		# self.E[1:-1,-1] = -self.state[1:-1,0]*(self.state[2:,-1]+self.state[:-2,-1]+self.state[1:-1,-2]+self.state[1:-1,0])
		# #top
		# self.E[-1,1:-1] = -self.state[-1,1:-1]*(self.state[0,1:-1]+self.state[-2,1:-1]+self.state[-1,:-2]+self.state[-1,2:])
		# #bottom
		# self.E[0,1:-1] = -self.state[0,1:-1]*(self.state[1,1:-1]+self.state[-1,1:-1]+self.state[0,:-2]+self.state[0,2:])

		# #bottom left
		# self.E[0,0] = -self.state[0,0]*(self.state[1,0]+self.state[-1,0]+self.state[0,-1]+self.state[0,1])
		# #bottom right
		# self.E[0,-1] = -self.state[0,-1]*(self.state[1,-1]+self.state[-1,-1]+self.state[0,-2]+self.state[0,0])
		# #top left
		# self.E[-1,0] = -self.state[-1,0]*(self.state[0,0]+self.state[-2,0]+self.state[-1,-1]+self.state[-1,1])
		# #top right
		# self.E[-1,-1] = -self.state[-1,-1]*(self.state[0,-1]+self.state[-2,-1]+self.state[-1,-2]+self.state[-1,0])

        return self.E


    def Choose_Spin(self):
        self.x = randint(0, self.L)
        self.y = randint(0, self.L)
        self.s = self.state[self.x, self.y]
        self.LE = self.E[self.x, self.y]
        return self.s , self.LE

    def Flip(self):
        if self.s == 1:
            self.state[self.x,self.y] = -1
        else:
            self.state[self.x,self.y] = 1
        self.fs = self.state[self.x,self.y]
        return self.fs

    def New_Energy(self):
        if self.fs == 1:
            self.NE = -self.fs*self.LE 
        else:
            self.NE = self.fs*self.LE
        return self.NE


    def Change_Energy(self):
        self.dE = self.NE - self.LE
        return self.dE


    def Decision(self):
        if self.dE < 0: 
            self.state[self.x,self.y] = self.fs
        elif uniform(0,1)< np.exp(-self.dE/self.T): 
            self.state[self.x,self.y] = self.fs
        else:
            self.state[self.x,self.y] = self.s 
        return self.state[self.x,self.y]

    def calcEnergy(self):
        energy = 0
        for i in range(self.L):
            for j in range(self.L):
                S = self.state[i,j]
                nb = self.state[(i+1)%self.L, j] + self.state[i,(j+1)%self.L] + self.state[(i-1)%self.L, j] + self.state[i,(j-1)%self.L]
                energy += -nb*S
                self.energy = energy/4.
        return self.energy

    def Magnetization(self):
        return np.sum(self.state)

    def MC_Step(self, step=10000):
        self.step = step 
        self.Energy1 = 0
        self.Energy2 = 0
        for _ in range(step):
            self.Hamiltonian()
            self.Choose_Spin()
            self.Flip()
            self.New_Energy()
            self.Change_Energy()
            self.Decision()
            self.Energy1 += self.calcEnergy()
            self.Energy2 += self.calcEnergy()**2
		
        return self.Magnetization(), self.Energy1, self.Energy2

if __name__ == "__main__":

    import os 
	
    try:
        os.mkdir("Data")
    except OSError:
        pass

	
    for temp in np.linspace(0.1, 4, 64):
        dt = (4 - 0.1) / 64
        y = Ising(20, temp)
        y.Initialize()
        Mag_Data = []
        cv_Data = []
        Energy1 = 0
        Energy2 = 0
        for i in range(0, 50, 2):	
            Mag, E1, E2= y.MC_Step()
            step = 10000
            cv = (E2 / step - E1*E1 / step**2) * dt**2 
            cv_Data.append(cv)
            Mag_Data.append(Mag)
			
	
		#Below command will save the data from the simulation 
        title_Mag = "Magnetization on running the simulation \n for x times under similar condition"
        title_cv = "Specific Heat on running the simulation \n for x times under similar condition"
        np.savetxt("Data/Mag/Temp_%s.csv" %temp , Mag_Data, fmt="%i", header = title_Mag)
        np.savetxt("Data/cv/Temp_%s.csv" %temp , cv_Data, fmt="%i", header = title_cv)
		
