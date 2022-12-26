# -*- coding: utf-8 -*-
"""
Created on Tue Oct 17 11:46:15 2017

@author: river603
"""

import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np
import csv
import os
import re
import glob

# path_list = 'E:/work/paper/ESPL_NumericalMorpho/RCEM/Run4/'
path_list = 'D:/HLL/kataoka/'
filename  = 'tt299'
#filename  = '2400min.txt'
#filename  = 'Run4_c6.txt'

#inputfiles = glob.glob(path_list+"*.txt")
inputfiles = glob.glob("*.txt")

#nx = 1000
#ny = 20

#nx = 180
#ny = 150

for inputfile in inputfiles:
    
#    filename = os.path.splitext(os.path.basename(inputfile))[0]
#    foldername = os.path.dirname(inputfile)
#    figname = foldername+"/fig/"+filename+".jpg"
    figname = "fig/"+inputfile+".jpg"

    with open(inputfile, 'r') as ifile:
        line = ifile.readline()
        head = line.replace('\n', '').split()
        
        tt = float(head[0])
        nx = int(head[1])
        ny = int(head[2])
        
        print(tt,nx,ny)
        
        lines = ifile.readlines()
            
        x  = []
        y  = []
        z  = []
        
        for line in lines:
           value = line.replace('\n', '').split()
           x.append(float(value[0]))
           y.append(float(value[1]))
           z.append(float(value[2]))
           
        
        X = [[0 for i in range(nx+1)] for j in range(ny+1)]
        Y = [[0 for i in range(nx+1)] for j in range(ny+1)]
        Z = [[0 for i in range(nx+1)] for j in range(ny+1)]
        
        ij = 0
    #    for i in range(nx+1):
    #        for j in range(ny+1):
        for j in range(ny+1):
            for i in range(nx+1):
                X[j][i] = x[ij]
                Y[j][i] = y[ij]-0.225
                Z[j][i] = z[ij]
                
                ij = ij+1
                
        
        plt.figure(figsize=(10,1.5))
        plt.axes().set_aspect('equal')
        plt.xlabel('X(m)')
        plt.ylabel('Y(m)')
        plt.subplots_adjust(bottom=0.3)
    #    plt.axis([30,45,0,0.9])
    #    plt.xticks( np.arange(30,50,5) )
        plt.axis([25,30,-0.25,0.25])
        plt.xticks( np.arange(25,31,1) )
        plt.yticks( np.arange(-0.25,0.3,0.25) )   
        levels = np.linspace(-0.02,0.01,51)
        labels = np.linspace(-0.02,0.01,6)
        
        plt.title("Time= {0:.0f} min".format(tt/60)) 
            
        plt.contourf(X, Y, Z, levels, cmap=cm.rainbow, extend="both")
        plt.colorbar(ticks=labels, label='Bed evolution (m)', orientation='horizontal', pad=0.1)
        plt.savefig(figname, dpi=600)
        plt.close()
            
    ifile.close()
        
        