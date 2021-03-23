#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 20 10:25:19 2020

@author: karthik & Pavi
"""


from CellType import CellType
import random as rand
import numpy as np
from tabulate import tabulate
import os
import sys
import matplotlib.pyplot as plt
from matplotlib import colors
import pandas as pd



class Simulate():
    def __init__(self, T, dT, XYMatrix, log = False):
        self.T = T
        self.dT = dT
        self.XYMatrix = XYMatrix
        self.XYFrames = []
        self.StoCRate = 0.
        self.CtoToxinRate = 0.008
        self.PlasmidPerCDeath = 20
        self.CtoDeathRate = 0.
        self.Cg1KillRate = 0.006
        self.Cg2KillRate = 0.002
        self.CtoCgRate = 0.
        self.CtoCg1Rate = 0.
        self.Cg1toCg2Rate = 0.
        self.Cg2toCg3Rate = 0.
        self.rCg = 0.0729
        self.KillProb = 5/10
        self.diffusivity = 5
        self.MeanVariation = 0.25
        self.K = 60000
        self.SCell = CellType("S Cell", 1, self.K, 0, 0)
        self.CCell = CellType("C Cell", 2, self.K, 0, 0)
        self.Cg1Cell = CellType("Cg1 Cell", 3, self.K, 0, 0)
        self.Cg2Cell = CellType("Cg2 Cell", 4, self.K, 0, 0)
        self.Cg3Cell = CellType("Cg3 Cell", 5, self.K, 0, 0)
        self.ToxinKills = 0 # USed for logging
        self.StoCConversions = 0 # Used for logging
        self.log = log
        
        self.tlist = [] #Used for logging
        self.SCelllist = []
        self.CCelllist = []
        self.Cg1Celllist = []
        self.Cg2Celllist = []
        self.Cg3Celllist = []
    
    def GenerateNo(self, Matrix, N, ID):
        n = 0
        while(n < N):
            Randomspot = np.random.choice(np.arange(int(0.1*len(Matrix)), len(Matrix)-int(0.1*len(Matrix))),2)
            if(Matrix[Randomspot[0], Randomspot[1]] == 0):
                Matrix[Randomspot[0], Randomspot[1]] = ID
                n+=1
            else:
                continue
        return Matrix
    
    def LogisticGrowth(self, NoOfCells, NoOfSCells, NoOfCCells, NoOfCg1Cells, NoOfCg2Cells, NoOfCg3Cells, r, K):
        N = NoOfCg1Cells + NoOfCg2Cells + NoOfCg3Cells + NoOfSCells + NoOfCCells
        return NoOfCells + NoOfCells*r*(1 - N/K)
    
    
    def Calculate_dN(self, Nt, NtMinus):
        return Nt - NtMinus
    

    def Transmutation(self, NoOfTransmutation, Matrix, DecayID, ID):
        NoOfConversions = NoOfTransmutation
        X,Y = np.where(Matrix==ID)
        if(len(X) > 0 and NoOfConversions > 0):
            index = np.random.choice(np.arange(0,len(X)), NoOfConversions)
            for i in index:
                Matrix[X[i]][Y[i]] = DecayID
        return Matrix



    def CheckNeighbours(self, Matrix, SeedX, SeedY, ID):
        choiceList = [(1,0), (0,1), (-1,0), (0,-1), (1,1), (-1,1), (-1,-1), (1,-1)]
        neighbourCoord = [-99,-99]
        NotFound = True
        while(NotFound and len(choiceList) > 0):
            randomNeighbourList = rand.choices(choiceList, k=1)
            randomNeighbour = randomNeighbourList[0]
            if(SeedX+randomNeighbour[0] >= len(Matrix) or SeedX+randomNeighbour[0] < 0 or SeedY+randomNeighbour[1] >= len(Matrix[0]) or SeedY+randomNeighbour[1] < 0):
                break
            elif(Matrix[SeedX+randomNeighbour[0],SeedY+randomNeighbour[1]] == 0):
                neighbourCoord = [SeedX+randomNeighbour[0], SeedY+randomNeighbour[1]]
                NotFound = False 
            elif(Matrix[SeedX+randomNeighbour[0],SeedY+randomNeighbour[1]] < ID and Matrix[SeedX+randomNeighbour[0],SeedY+randomNeighbour[1]] > 0):
                neighbourCoord = [SeedX+randomNeighbour[0], SeedY+randomNeighbour[1]]
                
                NotFound = False
                n_id = Matrix[SeedX+randomNeighbour[0],SeedY+randomNeighbour[1]]
                # print(n_id, neighbourCoord, self.tlist[-1]+1)
                Matrix = self.PutPoints(Matrix, 1, n_id )
                # print(n_id, neighbourCoord, self.tlist[-1]+1)
            else:
                choiceList.remove(randomNeighbour)

        return neighbourCoord, Matrix


    def DelElements(self, Matrix, ID, dN):
        X,Y = np.where(Matrix==ID)
        if(len(X) > 0):
            index = np.random.choice(np.arange(len(X)), dN)
            Matrix[X[index], Y[index]] = 0
        return Matrix    
    
    def PutPoints(self, Matrix, dN, ID):
        n = 0
        X,Y = np.where(Matrix == ID)
        X = list(X)
        Y = list(Y)
        while(n < dN and len(X) > 0):
            index = rand.randint(0, len(X)-1)
            SeedX = X[index]
            SeedY = Y[index]
            neighbourCoord, Matrix = self.CheckNeighbours(Matrix, SeedX, SeedY, ID)
            if(neighbourCoord[0] < 0):
                #print("Point Stuck in between all Cells. Calculating new centers for splitting ")
                X.pop(index)
                Y.pop(index)
                if(len(X)==0):
                    print("All the points of ID : {} has no free sites to split dN = {}; n = {}".format(ID, dN, n))

            else:
                Matrix[neighbourCoord[0], neighbourCoord[1]] = ID
                n+=1
                X.append(neighbourCoord[0])
                Y.append(neighbourCoord[1])
    
        return Matrix

    
        
    
    def GhostFunc(self, a, size):
        if(a<0):
            return 0
        elif(a>=size):
            return size - 1
        else:
            return a
    
    def CoordinateReturn(self, size, x, y, diffusivity):
        leftend = self.GhostFunc(x-diffusivity, size)
        rightend = self.GhostFunc(x+diffusivity, size)
        bottom = self.GhostFunc(y-diffusivity, size)
        top = self.GhostFunc(y+diffusivity, size)
        return leftend, rightend, bottom, top
    
    
    def KillSCells(self, Matrix, X, Y, ID):
        diffusivity = self.diffusivity       
        for x, y in zip(X,Y):
            leftend, rightend, bottom, top = self.CoordinateReturn(len(Matrix), x, y, diffusivity)
            Sindex_X, SindexY = np.where(Matrix[leftend:rightend+1,bottom:top+1]==1)
            
            for i,j in zip(leftend+Sindex_X, bottom+SindexY):
                if(rand.random() <= (self.KillProb)**(max(abs(i-x), abs(j-y)))):
                    Matrix[i,j] = 0
                    self.ToxinKills+=1

        return Matrix
    
    def ToxinDeaths(self, NoOfKills, Matrix, ID, StoCCellLambda = 0.):
        X,Y = np.where(Matrix==ID)
        if(len(X) <= 0 or NoOfKills <= 0):
            return Matrix
        else:
            index = np.random.choice(np.arange(0,len(X)), NoOfKills)
            Matrix[X[index], Y[index]] = 0 # setting the Toxin producer to 0 i.e it dies
            if(ID == 2 and StoCCellLambda > 0):
                for x, y in zip(X[index],Y[index]): # loop through selected ID=3 locations from index
                    leftend, rightend, bottom, top = self.CoordinateReturn(len(Matrix), x, y, self.diffusivity)
                    temp_NoOfSCells = np.count_nonzero(Matrix[leftend:rightend+1, bottom:top+1] == 1) # counting the number of S cells ID==1 around each dead C cells
                    NoOfConversions = np.sum(np.random.binomial(temp_NoOfSCells, StoCCellLambda, rand.randint(self.PlasmidPerCDeath - 5, self.PlasmidPerCDeath))) # Need to improve upon this
                    #if(temp_NoOfSCells >= 110): print("No of SCells found is : {} \n No of Conversions is : {}".format(temp_NoOfSCells, NoOfConversions))
                    #print("No Of Conversions from S to C with plasmid is : {}".format(NoOfConversions))
                    if(NoOfConversions > 0 and temp_NoOfSCells > 0):
                        SindexX, SindexY = np.where(Matrix[leftend:rightend+1, bottom:top+1]==1)
                        placepos = np.random.choice(np.arange(0, len(SindexX)), NoOfConversions)
                        #print(Matrix[leftend+SindexX[placepos],bottom+SindexY[placepos]])
                        Matrix[leftend+SindexX[placepos],bottom+SindexY[placepos]] = 2
                        self.StoCConversions+=NoOfConversions
                        
            
            Matrix = self.KillSCells(Matrix, X[index], Y[index], ID)
               # Matrix = self.KillSCells(Matrix, X[index], Y[index], 2, 0.005)
            return Matrix


    def logger(self, name, delT, dNS, dNC, dNCg1, dNCg2, dNCg3, N_CtoToxin, N_CtoDeath, N_CtoCg, cg1, cg2, cg3, NoOfToxinKills, N_StoC, N_Cg1Kills, N_Cg2Kills):
        
        file = open(name + ".csv", "a")
        if(delT == 1):
            file.write("TimeStep,SCells,CCells,Cg1Cells,Cg2Cells,Cg3Cells,dNS,dNC,dNCg1,dNCg2,dNCg3,NoOfCtoToxin,NoOfCtoDeath,NoOfCtoCg(NotUsed),NoOfCtoCg1,NoOfCg1toCg2,NoOfCg2toCg3,NoOfToxinKills,NoOfStoC,NoOfCg1Kills,NoOfCg2Kills\n")
        file.write(str(delT) + "," + str(self.SCell.NtMinus) + "," + str(self.CCell.NtMinus) + "," + str(self.Cg1Cell.NtMinus) + "," + str(self.Cg2Cell.NtMinus) + "," + str(self.Cg3Cell.NtMinus) + "," + str(dNS) + "," + str(dNC) 
                   + "," + str(dNCg1) + "," + str(dNCg2) + "," + str(dNCg3)
                    + "," + str(N_CtoToxin) + "," + str(N_CtoDeath) + "," + str(N_CtoCg) + "," + str(cg1) + "," + str(cg2) + "," + str(cg3) + "," + str(NoOfToxinKills) 
                    + "," + str(N_StoC) + "," + str(N_Cg1Kills) + "," + str(N_Cg2Kills) + "\n")
        
        file.close()
        
    
    
    def PrintInfo(self, Name):
        with open(Name+"Info.txt","w") as f:
            
            f.write("############## Starting the Simulation #################### \n \n")
            f.write("Parameters that are used for this simulation are: \n")
            data = [["Total Time ", self.T, 0], ["Time Steps ", self.dT, 0], 
                    ["Grid Size", np.shape(self.XYMatrix), 0], 
                    ["No Of Intial S Cells", self.SCell.N0, 0],
                    ["No Of Initial C Cells", self.CCell.N0, 0], 
                    ["Growth Rate for S Cells", self.SCell.r, 0],
                    ["Growth Rate for C Cells", self.CCell.r, 0],
                    ["Lambda S to C", self.StoCRate, self.MeanVariation*self.StoCRate],
                    ["Lambda C to ToxinRate", self.CtoToxinRate, self.MeanVariation*self.CtoToxinRate],
                    ["C to Death Rate", self.CtoDeathRate, 0],
                    ["C to Cg Rate", self.CtoCgRate, self.MeanVariation*self.CtoCgRate],
                    ["Growth rate for Cg", self.Cg1Cell.r ,0],
                    ["Cg1 to Toxin rate", self.Cg1KillRate, self.MeanVariation*self.Cg1KillRate],
                    ["Cg2 to Toxin rate", self.Cg2KillRate, self.MeanVariation*self.Cg2KillRate],
                    ["Kill probability", self.KillProb, 0],
                    ["Kill Radius (both sides)", self.diffusivity*2,0]]
            f.write(tabulate(data, headers=["Parameter Names", "Values", "Variations"]))
    
    
    def ReturnWholeNumber(self, dN, decimal):
        if(decimal>=1):
            return int(dN)+int(decimal), decimal-int(decimal)
        else:
            return int(dN), decimal + dN - int(dN)
    
    def saveframes(self, XYFrames, t):
        if(t%5!=0):
            return 0
        if(os.path.exists("Images")==False):
            os.mkdir("Images")
        
        foldername = "Images"
        
        fig, ax = plt.subplots(1,2,figsize=(18,6), gridspec_kw={'width_ratios': [2, 1]})
        plt.ion()
        names = ["S Cells", "C Cells", "Cg1 Cells", "Cg2 Cells", "Cg3 Cells"]
        colors1 = ["r-", "y-", "g-", "b-", "c-", "m-"]
        cmap = colors.ListedColormap(['k', 'red', "yellow", "green", "blue", "cyan"])
        
        
        ax[0].plot(self.tlist, self.SCelllist, colors1[0], label = "{}".format(names[0]))
        ax[0].plot(self.tlist, self.CCelllist, colors1[1], label = "{}".format(names[1]))
        ax[0].plot(self.tlist, self.Cg1Celllist, colors1[2], label = "{}".format(names[2]))
        ax[0].plot(self.tlist, self.Cg2Celllist, colors1[3], label = "{}".format(names[3]))
        ax[0].plot(self.tlist, self.Cg3Celllist, colors1[4], label = "{}".format(names[4]))
        
        ax[0].set_xlim(0, self.T)
        ax[0].set_ylim(0, self.K + 750)
        ax[0].set_xlabel("Time in [arb units]", fontsize = 15)
        ax[0].set_ylabel("Population", fontsize = 15)
        ax[0].legend(fontsize = 20)
        mat = ax[1].imshow(XYFrames, origin="lower", cmap=cmap, vmin=0, vmax=5)
        cbar = plt.colorbar(mat, ax=ax[1], ticks = [0,1,2,3,4,5], fraction=0.046, pad=0.04)
        cbar.ax.set_yticklabels([" ", "S", "C", "Cg1", "Cg2", "Cg3"])
        fig.tight_layout()
        plt.ioff()
        plt.savefig(os.path.join(os.getcwd(),foldername,"frame{}.png".format(t)),dpi=50)
        
        plt.close(fig)    
    
    def write_csv(self, name):
        Cells = {"SCells" : self.SCelllist, "CCells" : self.CCelllist, "Cg1Cells" : self.Cg1Celllist, "Cg2Cells" : self.Cg2Celllist, "Cg3Cells" : self.Cg3Celllist}
        pd.DataFrame(data = Cells).to_csv(name + ".csv", index = False)
            
    
    
    def StartSimulation(self):

        print("############## Starting the Simulation ##################### \n \n ")
        print('''The Parameters used for running the simulation are: \n
        Total Time = {} \n
        Time steps = {} \n
        Matrix shape is = {} \n '''.format(self.T, self.dT, np.shape(self.XYMatrix)))

        XYMatrix = self.XYMatrix # Setting the Grid all zeros to start with
        loggername = "Logger"
        if(self.log and os.path.isfile("Logger.csv")):
            x = input("Removing the old log file is that ok ? yes/no : ")
            if(x=="yes"):
                os.remove("Logger.csv")
            else:
                sys.exit(0)
        SPopCounter = 0.0
        CPopCounter = 0.0
#        StoCCounter = 0.0 # There is No S to C at every step
        CtoToxinCounter = 0.0
        CDeathCounter = 0.0
        Cg1PopCounter = 0.0
        Cg2PopCounter = 0.0
        Cg3PopCounter = 0.0
        CtoToxinCounter = 0.0
        CtoCgCounter = 0.0
        Cg1KillCounter = 0.0
        Cg2KillCounter = 0.0

        # Defining the 0th time parameters. Starting with SCell
        
        self.SCell.Nt = self.SCell.N0 # Setting Nt to N0. No logic here though
        XYMatrix = self.GenerateNo(XYMatrix, self.SCell.Nt, self.SCell.ID) # Generate Nt==N0 cells in the Grid
        self.SCell.NtMinus = self.SCell.Nt # Setting NtMinus to current Nt
        
        
        self.CCell.Nt = self.CCell.N0 # Setting Nt to N0. No logic here though
        XYMatrix = self.GenerateNo(XYMatrix, self.CCell.Nt, self.CCell.ID) # Generate Nt==N0 cells in the Grid
        self.CCell.NtMinus = self.CCell.Nt # Setting NtMinus to Current Nt
        #self.XYFrames.append(XYMatrix.copy()) #Append the XYGrid at t0 into the frame list
        
        self.tlist.append(0)
        self.SCelllist.append(self.SCell.NtMinus)
        self.CCelllist.append(self.CCell.NtMinus)
        self.Cg1Celllist.append(self.Cg1Cell.NtMinus)
        self.Cg2Celllist.append(self.Cg2Cell.NtMinus)
        self.Cg3Celllist.append(self.Cg3Cell.NtMinus)
        
        
        if(self.log): 
            self.saveframes(XYMatrix, 0)
        
        for t in range(1, int(self.T/self.dT)):
            
            
            if(t%10==0): print("\n Running the time step {}".format(t))
            
            
    
############## Calculate the SCell Dynamics ###################################
            
            
            
            # Logistic Growth of S
            
            
            self.SCell.Nt = self.LogisticGrowth(self.SCell.NtMinus, self.SCell.NtMinus, 
                                                self.CCell.NtMinus, self.Cg1Cell.NtMinus, 
                                                self.Cg2Cell.NtMinus, self.Cg3Cell.NtMinus,
                                                self.SCell.r, self.SCell.K) # Now grow the Population
            
            dN_S = self.Calculate_dN(self.SCell.Nt, self.SCell.NtMinus) #Calculate the difference between Nt and NtMinus
            dN_SInt, SPopCounter = self.ReturnWholeNumber(dN_S, SPopCounter)
            
            
            # Transmutation from S to C
            # This is commented since there is no conversions from S to C at a constant rate
            #tempStoCRate = abs(rand.gauss(self.StoCRate, self.StoCRate*self.MeanVariation))
#            N_StoC = tempStoCRate*self.SCell.NtMinus
#            N_StoC, StoCCounter = self.ReturnWholeNumber(N_StoC, StoCCounter)
#            print(N_StoC)

###################### Cg Cell dynamics ####################
            
            # Toxin producing Cg1 and Cg2 Cells
            

            #tempCg1KillRate = abs(rand.gauss(self.Cg1KillRate, self.Cg1KillRate*self.MeanVariation))
            #N_Cg1Kill = tempCg1KillRate*self.Cg1Cell.NtMinus
            #N_Cg1Kill, Cg1KillCounter = self.ReturnWholeNumber(N_Cg1Kill, Cg1KillCounter)
            N_Cg1Kill = np.random.binomial(self.Cg1Cell.NtMinus, self.Cg1KillRate)
            N_Cg2Kill = np.random.binomial(self.Cg2Cell.NtMinus, self.Cg2KillRate)
            
            if(N_Cg1Kill >= 1):
                XYMatrix = self.ToxinDeaths(N_Cg1Kill, XYMatrix, self.Cg1Cell.ID)
            
            #tempCg2KillRate = abs(rand.gauss(self.Cg2KillRate, self.Cg2KillRate*self.MeanVariation))
            #N_Cg2Kill = tempCg2KillRate*self.Cg2Cell.NtMinus
            #N_Cg2Kill, Cg2KillCounter = self.ReturnWholeNumber(N_Cg2Kill, Cg2KillCounter)
            
            if(N_Cg2Kill >= 1):
                XYMatrix = self.ToxinDeaths(N_Cg2Kill, XYMatrix, self.Cg2Cell.ID)
            
            self.Cg1Cell.Nt = self.LogisticGrowth(self.Cg1Cell.NtMinus, self.SCell.NtMinus, 
                                                self.CCell.NtMinus, self.Cg1Cell.NtMinus, 
                                                self.Cg2Cell.NtMinus, self.Cg3Cell.NtMinus,
                                                self.Cg1Cell.r, self.Cg1Cell.K) # Now grow the Population
            
            dN_Cg1 = self.Calculate_dN(self.Cg1Cell.Nt, self.Cg1Cell.NtMinus)
            dN_Cg1Int, Cg1PopCounter = self.ReturnWholeNumber(dN_Cg1, Cg1PopCounter)
            
            # Grow Cg2 Cells
            
            self.Cg2Cell.Nt = self.LogisticGrowth(self.Cg2Cell.NtMinus, self.SCell.NtMinus, 
                                                self.CCell.NtMinus, self.Cg1Cell.NtMinus, 
                                                self.Cg2Cell.NtMinus, self.Cg3Cell.NtMinus,
                                                self.Cg2Cell.r, self.Cg2Cell.K) # Now grow the Population
            
            dN_Cg2 = self.Calculate_dN(self.Cg2Cell.Nt, self.Cg2Cell.NtMinus)
            dN_Cg2Int, Cg2PopCounter = self.ReturnWholeNumber(dN_Cg2, Cg2PopCounter)
            
            
            # Grow Cg3 Cells

            self.Cg3Cell.Nt = self.LogisticGrowth(self.Cg3Cell.NtMinus, self.SCell.NtMinus, 
                                                self.CCell.NtMinus, self.Cg1Cell.NtMinus, 
                                                self.Cg2Cell.NtMinus, self.Cg3Cell.NtMinus,
                                                self.Cg3Cell.r, self.Cg3Cell.K) # Now grow the Population
            
            dN_Cg3 = self.Calculate_dN(self.Cg3Cell.Nt, self.Cg3Cell.NtMinus)
            dN_Cg3Int, Cg3PopCounter = self.ReturnWholeNumber(dN_Cg3, Cg3PopCounter)            
            
            
            #tempCtoCgRate = abs(rand.gauss(self.CtoCgRate, self.CtoCgRate*self.MeanVariation))
            #N_CtoCg = tempCtoCgRate*self.CCell.NtMinus
            #N_CtoCg, CtoCgCounter = self.ReturnWholeNumber(N_CtoCg, CtoCgCounter)
            N_CtoCg = np.random.binomial(self.CCell.NtMinus, self.CtoCgRate)
            TotalNoOfCgTypes = np.random.choice([4,5,6], N_CtoCg)

            # Now randomly decay the Cg to Cg1 or Cg2 or Cg3

            #Cg1Cells = np.count_nonzero(TotalNoOfCgTypes==4) # No of Cg1Cells that are formed from C

            #Cg2Cells = np.count_nonzero(TotalNoOfCgTypes==5) # No of Cg2 Cells that are formed from C
            
            #Cg3Cells = np.count_nonzero(TotalNoOfCgTypes==6) # No of Cg3 Cells that are formed from C
                    
            Cg1Cells = np.random.binomial(self.CCell.NtMinus, self.CtoCg1Rate) # No of Cg1Cells that are formed from C

            Cg2Cells = np.random.binomial(self.Cg1Cell.NtMinus, self.Cg1toCg2Rate) # No of Cg2 Cells that are formed from C
            
            Cg3Cells = np.random.binomial(self.Cg2Cell.NtMinus, self.Cg2toCg3Rate) # No of Cg3 Cells that are formed from C
                    
                     
            
            
############## Calculate the CCell Dynamics ###################################            
            
            # Some of the CCells converts to Toxin and kills
            # Transmutation from C to Toxin
            
            #tempCtoToxinRate = abs(rand.gauss(self.CtoToxinRate, self.CtoToxinRate*self.MeanVariation))
            #N_CtoToxin = tempCtoToxinRate*self.CCell.NtMinus
            #N_CtoToxin, CtoToxinCounter = self.ReturnWholeNumber(N_CtoToxin, CtoToxinCounter)
            N_CtoToxin = np.random.binomial(self.CCell.NtMinus, self.CtoToxinRate)
            XYMatrix = self.ToxinDeaths(N_CtoToxin, XYMatrix, self.CCell.ID, self.StoCRate)
            
            # Death of C Cells
            
            #tempCDeathRate = abs(rand.gauss(self.CtoDeathRate, self.CtoDeathRate*self.MeanVariation))
            #N_CDeaths = tempCDeathRate*self.CCell.NtMinus
            #N_CDeaths, CDeathCounter = self.ReturnWholeNumber(N_CDeaths, CDeathCounter)
            N_CDeaths = np.random.binomial(self.CCell.NtMinus, self.CtoDeathRate)
            
            
            if(N_CDeaths >= 1): # Delete the C Cells if it is greater than 1
                XYMatrix = self.DelElements(XYMatrix, self.CCell.ID, N_CDeaths)            
            
            # Logistic Growth of C
            
            self.CCell.Nt = self.LogisticGrowth(self.CCell.NtMinus, self.SCell.NtMinus, 
                                                self.CCell.NtMinus, self.Cg1Cell.NtMinus, 
                                                self.Cg2Cell.NtMinus, self.Cg3Cell.NtMinus,
                                                self.CCell.r, self.CCell.K) # Now grow the Population
            
            
            dN_C = self.Calculate_dN(self.CCell.Nt, self.CCell.NtMinus) #Calculate the difference between Nt and NtMinus
            dN_CInt, CPopCounter = self.ReturnWholeNumber(dN_C, CPopCounter)        
            

######################## Filling in the Matrix for conversions and logistic Growths ##################
            
#            XYMatrix = self.Transmutation(N_StoC, XYMatrix, DecayID = self.CCell.ID, ID = self.SCell.ID)
            if(dN_SInt>0):
                XYMatrix = self.PutPoints(XYMatrix, dN_SInt, self.SCell.ID) # add points if dN is positive
            else:
                XYMatrix = self.DelElements(XYMatrix, self.SCell.ID, -1*dN_SInt) # remove points if dN is negative
            
   
            if(dN_CInt>0):
                XYMatrix = self.PutPoints(XYMatrix, dN_CInt, self.CCell.ID) # add points if dN is positive
            else:
                XYMatrix = self.DelElements(XYMatrix, self.CCell.ID, -1*dN_CInt) # remove points if dN is negative
             
            
            
            XYMatrix = self.Transmutation(Cg1Cells, XYMatrix, self.Cg1Cell.ID, self.CCell.ID)
            if(dN_Cg1Int > 0):
                XYMatrix = self.PutPoints(XYMatrix, dN_Cg1Int, self.Cg1Cell.ID) # add points if dN is positive
            else:
                XYMatrix = self.DelElements(XYMatrix, self.Cg1Cell.ID, -1*dN_Cg1Int) # remove points if dN is negative
            
            
            XYMatrix = self.Transmutation(Cg2Cells, XYMatrix, self.Cg2Cell.ID, self.Cg1Cell.ID)
            if(dN_Cg2Int > 0):
                XYMatrix = self.PutPoints(XYMatrix, dN_Cg2Int, self.Cg2Cell.ID) # add points if dN is positive
            else:
                XYMatrix = self.DelElements(XYMatrix, self.Cg2Cell.ID, -1*dN_Cg2Int) # remove points if dN is negative
            
            
            XYMatrix = self.Transmutation(Cg3Cells, XYMatrix, self.Cg3Cell.ID, self.Cg2Cell.ID)
            if(dN_Cg3Int > 0):
                XYMatrix = self.PutPoints(XYMatrix, dN_Cg3Int, self.Cg3Cell.ID) # add points if dN is positive
            else:
                XYMatrix = self.DelElements(XYMatrix, self.Cg3Cell.ID, -1*dN_Cg3Int) # remove points if dN is negative
                        
    
            
###################### Calculate all the NtMinus for S, C, Cg1, Cg2, Cg3 cells #######################                                    
            
            
            
            self.SCell.NtMinus = np.count_nonzero(XYMatrix==self.SCell.ID) + dN_S - dN_SInt
            self.CCell.NtMinus = np.count_nonzero(XYMatrix==self.CCell.ID) + dN_C - dN_CInt
            self.Cg1Cell.NtMinus = np.count_nonzero(XYMatrix==self.Cg1Cell.ID) + dN_Cg1 - dN_Cg1Int
            self.Cg2Cell.NtMinus = np.count_nonzero(XYMatrix==self.Cg2Cell.ID) + dN_Cg2 - dN_Cg2Int
            self.Cg3Cell.NtMinus = np.count_nonzero(XYMatrix==self.Cg3Cell.ID) + dN_Cg3 - dN_Cg3Int
            
            
            self.tlist.append(t)
            self.SCelllist.append(int(self.SCell.NtMinus))
            self.CCelllist.append(int(self.CCell.NtMinus))
            self.Cg1Celllist.append(int(self.Cg1Cell.NtMinus))
            self.Cg2Celllist.append(int(self.Cg2Cell.NtMinus))
            self.Cg3Celllist.append(int(self.Cg3Cell.NtMinus))
            
            #self.XYFrames.append(XYMatrix.copy())
            
            if(self.log):
                self.logger(loggername, t, dN_S, dN_C, dN_Cg1, dN_Cg2, dN_Cg3, N_CtoToxin, N_CDeaths, N_CtoCg, Cg1Cells, Cg2Cells, Cg3Cells, self.ToxinKills, self.StoCConversions, N_Cg1Kill, N_Cg2Kill)
                self.saveframes(XYMatrix, t)
                        
            self.ToxinKills = 0 # resetting the logger variable
            self.StoCConversions = 0 #resetting the logger variable
    
    
    
