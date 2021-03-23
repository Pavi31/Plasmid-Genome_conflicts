#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Nov 21 13:46:24 2020

@author: karthik
"""


import matplotlib.pyplot as plt
import numpy as np
import matplotlib
import matplotlib.animation as animation
import Simulation
from IPython.display import HTML
# import plotly.graph_objects as go
# from plotly.offline import plot
from matplotlib import colors

from matplotlib.animation import FuncAnimation, PillowWriter 


def ShowGrid(Matrix):
    plt.figure(figsize=(10,10))
    plt.imshow(Matrix, origin = "lower", cmap = "plasma")
    plt.colorbar()
    plt.show()
    
def PlotTimeLine(ListXYFrames, rows, columns):
    fig, ax = plt.subplots(rows, columns, figsize=(12,12), sharex = True, sharey = True)
    intervals = int(len(ListXYFrames)/(rows*columns))
    counts = 0
    cmap = colors.ListedColormap(['k', 'red', "yellow", "green", "blue", "cyan"])
    for i in range(rows):
        for j in range(columns):
            im = ax[i,j].imshow(ListXYFrames[counts*intervals], origin="lower", cmap=cmap, vmin=0, vmax=5)
            ax[i,j].get_xaxis().set_visible(0)
            ax[i,j].get_yaxis().set_visible(0)
            ax[i,j].set_title("Time t = {}".format(counts*intervals))
            counts+=1
    fig.colorbar(im, ax=ax.ravel().tolist())
    plt.show()

def CellPop(SimulationObj):
    plt.figure(figsize = (12,6))
    for i in SimulationObj.ListCellType:
        i.PlotPop()
    plt.xlabel("Time in [arb units]", fontsize = 15)
    plt.ylabel("Population", fontsize = 15)
    plt.legend(fontsize = 20)
    plt.show()
    
def CellPopGrid(SimulationObj, NoOfID):
    names = ["S Cells", "C Cells", "Cg1 Cells", "Cg2 Cells", "Cg3 Cells"]
    colors = ["r-", "y-", "g-", "b-", "c-", "m-"]
    plt.figure(figsize = (12,6))
    for ids in range(1,NoOfID+1):
        temp_lst = [np.count_nonzero(i == ids) for i in SimulationObj.XYFrames]
        plt.plot(temp_lst, colors[ids-1], label = "{}".format(names[ids-1]))
    plt.xlabel("Time in [arb units]", fontsize = 15)
    plt.ylabel("Population", fontsize = 15)
    plt.legend(fontsize = 20)
    plt.show()
def Animation(ListXYGrid):
    fig, ax = plt.subplots()
    def animate(i):
        im = ax.imshow(ListXYGrid[i], origin="lower", cmap=plt.cm.get_cmap('jet', 7), vmin=0, vmax=6)
    ani = matplotlib.animation.FuncAnimation(fig, animate, len(ListXYGrid))
    HTML(ani.to_jshtml())
    ani.save('animation1.gif', writer='imagemagick', fps=100)

def plotall(XYFrames):
    plt.figure(figsize = (12,6))
    for i in range(1,4):
        plt.plot([np.count_nonzero(j==i) for j in XYFrames], label = "Type {}".format(i))
    plt.legend()
    plt.show()


def PlotTrialTimeLine(SimulationList, IDS):
    AvgPop = [[], [], [], [], [], []]
    AvgPopError = [[], [], [], [], [], []]
    colors = ["red", "green", "blue", "cyan", "yellow", "magenta"]
    names = ["S Cells", "C Cells", "COn Cells", "Cg1 Cells", "Cg2 Cells", "Cg3 Cells"]
    plt.figure(figsize = (12,6))
    for time in range(len(SimulationList[0].XYFrames)):
        tempPop = [[],[],[],[],[],[]]
        for n in SimulationList:
            for ids in range(1, IDS+1):
                tempPop[ids-1].append(np.count_nonzero(n.XYFrames[time] == ids))
        
        for temp in range(IDS):
            tempPop[temp] = np.array(tempPop[temp])
            AvgPop[temp].append(tempPop[temp].mean())
            AvgPopError[temp].append(tempPop[temp].std()/np.sqrt(len(SimulationList)))
    for i in range(len(AvgPop)):
        plt.plot(np.arange(0,len(AvgPop[0])), AvgPop[i], label = "{}".format(names[i]), linewidth = 2)
        plt.fill_between(np.arange(0,len(AvgPop[0])), np.array(AvgPop[i]) - np.array(AvgPopError[i]),
                         np.array(AvgPop[i]) + np.array(AvgPopError[i]), alpha=0.5, linewidth=0)
    plt.legend(fontsize = 15)
    plt.show()

def Animate(XYFrames):
    fig, ax = plt.subplots(figsize = (12, 12))
    ims = []
    def update(i):
        cmap = colors.ListedColormap(['k', 'red', "yellow", "green", "blue", "cyan", "m"])
        ax.imshow(XYFrames[i], origin = "lower", cmap=cmap, vmin=0, vmax=5)
        ims.append(XYFrames)
    ani = FuncAnimation(fig, ims, interval=50, blit=True,
                                repeat_delay=1000)
    plt.show()