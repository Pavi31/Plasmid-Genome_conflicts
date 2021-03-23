#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 20 10:32:35 2020

@author: karthik & Pavi
"""

import Simulation
import Lattice
from CellType import CellType
import numpy as np
import Visualisation as Vis
import pandas as pd
import os
#import VideoMaker


def WriteToCSV(XYFrames, name):
    Cells = {"SCells" : [], "CCells" : [], "Cg1Cells" : [], "Cg2Cells" :[], "Cg3Cells" : []}
    for i in range(len(XYFrames)):
        for ids, key in enumerate(Cells.keys()):
            Cells[key].append(np.count_nonzero(XYFrames[i] == ids+1))
    pd.DataFrame(data = Cells).to_csv(name + ".csv", index = False)




for i in range(1,11):
    FolderName = 'DraftIter' + str(i)
    
    os.mkdir(FolderName)
    os.chdir(FolderName)
    
    XYGrid = Lattice.Lattice(300, 300).CreateLattice()
    
    Model = Simulation.Simulate(10000, 1, XYGrid, log = True)
    
    print(np.ndim(Model.XYMatrix))
    Model.SCell = CellType("S Cell", ID = 1, K = 60000, r = 0.0729, N0 = 200)
    
    Model.CCell = CellType("C Cell", ID = 2, K = 60000, r = 0.0729, N0 = 2)
    
    Model.Cg1Cell = CellType("Cg1 Cell", ID = 3, K = 60000, r = 0.0729, N0 = 0)
    
    Model.Cg2Cell = CellType("Cg2 Cell", ID = 4, K = 60000, r = 0.0729, N0 = 0)
    
    Model.Cg3Cell = CellType("Cg3 Cell", ID = 5, K = 60000, r = 0.0729, N0 = 0)
    
    Model.StoCRate = 0.001 # sweeping this
    Model.CtoToxinRate = 0.008
    Model.PlasmidPerCDeath = 20 # currently vary between PlasmidPerCDeath - 5 to PlasmidPerCDeath
    Model.CtoDeathRate = 0.0001 # sweeping this
    Model.Cg1KillRate = 0.006
    Model.Cg2KillRate = 0.002
    Model.CtoCg1Rate = 0.000001 # sweeping this
    Model.Cg1toCg2Rate = 0.0001
    Model.Cg2toCg3Rate = 0.0001
    Model.rCg = 0.0729
    Model.KillProb = 5/10 # sweeping this
    Model.diffusivity = 5
    Model.MeanVariation = 0.0
    Model.K = 60000
    
    Model.PrintInfo("test")
    
    Model.StartSimulation()
    Model.write_csv(FolderName)
    os.chdir(os.getcwd()+"/../")
    if(Model.log):
        path = os.path.join(FolderName,"Images")
        VideoMaker.generate_video(path, i)
    del Model
