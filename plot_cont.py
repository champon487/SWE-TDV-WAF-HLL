# -*- coding: utf-8 -*-
"""
Created on Tue Oct 17 11:46:15 2017

@author: river603
"""

from numba import njit, jit, f8
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np
import csv
import os
import re
import glob
import shutil


def length_height(x,y,z,nnn):
    
    s  = x[0,:]
    wdz = z[0,:]
    dz = wdz-np.mean(wdz)
    
    window = 5 # 移動平均の範囲
    w = np.ones(window)/window

    dz2 = np.convolve(dz, w, mode='same')
        
    
    zmin = []
    zmax = []
    hb = []
    
    nx = len(s)-1
    
    for i in range(nx+1):
        zmin.append(min(z[:,i]))
        zmax.append(max(z[:,i]))
        hb.append(max(z[:,i])-min(z[:,i]))
#        print(s[i],min(z[:,i]),max(z[:,i]))
        
    ip = []
    
    for i in range(nx):
        if np.sign(dz2[i])<0. and np.sign(dz2[i+1])>0. and np.abs(zmax[i]-zmin[i])>0.001:
            ip.append(i)
                
    nn = len(ip)-1
    
    ipp = 0
    
#    if np.mod(nn,2)!=0:
#        ipp = 1
    
    ll = []
    hh = []
    
#    for n in range(0,nn-ipp,2):
    for n in range(0,nn-1):
        wl = s[ip[n+1]]-s[ip[n]]
#        wh = np.max(zmax[ip[n]:ip[n+1]])-np.min(zmin[ip[n]:ip[n+1]])
        wh = np.max(hb[ip[n]:ip[n+1]])
        
#        if 80>s[ip[n]]>40.:
        if 110>s[ip[n]]>80.:
        
            ll.append(wl)
            hh.append(wh)

    L = 0.
    H = 0.
    
#    nb = len(ll)
    
#    if nb>=1:
            
#        L = np.mean(ll)
#        H = np.mean(hh)
                
#    return L, H
#    plt.figure(figsize=(10,2.5),tight_layout=True)
#    plt.ylim([-0.05,0.025])
#    plt.plot(s,zmin)
#    plt.plot(s,zmax)
#    plt.plot(s,dz2)
    
#    for n in range(0,nn-1):
#        plt.plot(s[ip[n]],dz2[ip[n]],'o')
        
#    plt.savefig('bar'+str(nnn)+'.jpg', dpi=600)
#    plt.close()
    

    return ll, hh    
    
    
#path_list = 'D:/HLL/original/aaaa/akahori/f-sp-bed1upwind/'
#path_list = 'D:/HLL/original/aaaa/wata/f-van-1st-am075/'
path_list = 'D:/HLL/original/aaaa/wata/f-hll/'
#path_list = 'D:/HLL/original/aaaa/fine/'
#path_list = 'D:/HLL/original/watanabeS-50/superbee/coarse/2ndorder-fix2/'
#path_list = 'D:/HLL/original/akahori/superbee/corase-fix/'
filename  = 'tt299'
#filename  = '2400min.txt'
#filename  = 'Run4_c6.txt'

inputfiles = glob.glob(path_list+"*.txt")
#inputfiles = glob.glob("*.txt")

fig_path = path_list+"/fig/"
if os.path.exists(fig_path):
    shutil.rmtree(fig_path)

os.mkdir(fig_path)

#nx = 1000
#ny = 20

#nx = 180
#ny = 150

output = open(path_list+"T-LH.dat",'w')

output.write("Time(min)  Length(m)  Height(m) \n")

nn = 0


for inputfile in inputfiles[::2]:
    
    filename = os.path.splitext(os.path.basename(inputfile))[0]
    foldername = os.path.dirname(inputfile)
    figname = foldername+"/fig/"+filename+".jpg"
#    figname = "fig/"+inputfile+".jpg"

    with open(inputfile, 'r') as ifile:
        line = ifile.readline()
        head = line.replace('\n', '').split()
        
        tt = float(head[0])
        nx = int(head[1])
        ny = int(head[2])
                
        lines = ifile.readlines()
            
        x  = []
        y  = []
        z  = []
        
        for line in lines:
           value = line.replace('\n', '').split()
           x.append(float(value[0]))
           y.append(float(value[1]))
           z.append(float(value[2]))
           
        
        wX = [[0 for i in range(nx+1)] for j in range(ny+1)]
        wY = [[0 for i in range(nx+1)] for j in range(ny+1)]
        wZ = [[0 for i in range(nx+1)] for j in range(ny+1)]
        
        X = np.array(wX, dtype="f8")
        Y = np.array(wY, dtype="f8")
        Z = np.array(wZ, dtype="f8")
        
        y_offset = 0.45
        
        ij = 0
    #    for i in range(nx+1):
    #        for j in range(ny+1):
        for j in range(ny+1):
            for i in range(nx+1):
                X[j][i] = x[ij]
                Y[j][i] = y[ij]-y_offset
                Z[j][i] = z[ij]
                
                ij = ij+1
                
        L, H = length_height(X,Y,Z,nn)
        
        if len(L)>0:
            lave = np.mean(L)
            have = np.mean(H)
            lmax = np.max(L)
            hmax = np.max(H)
            lmin = np.min(L)
            hmin = np.min(H)
        else:
            lave = 0.
            have = 0.
            lmax = 0.
            hmax = 0.
            lmin = 0.
            hmin = 0.
        
        print(tt/60,lave,have)
        
        
#        aaa = str(tt/60.)+"\t"+str(L)+"\t"+str(H)+"\n"
        aaa = str(tt/60.)+"\t"+str(lave)+"\t"+str(have)+"\t"+str(lmax)+"\t"+str(hmax)+"\t"+str(lmin)+"\t"+str(hmin)+"\n"
        output.write(aaa)
                
        plt.figure(figsize=(10,2.5),tight_layout=True)
        plt.axes().set_aspect('equal')
        plt.xlabel('X(m)')
        plt.ylabel('Y(m)')
        plt.subplots_adjust(bottom=0.2)
        
        zmin = -0.05
        zmax =  0.025
#        xmax = 115.
#        xmin = 100.
        xmax = 105.
        xmin = 90.
#        xmax = 140.
#        xmin = 80.
        ymin = -y_offset
        ymax =  y_offset
        
    #    plt.axis([30,45,0,0.9])
    #    plt.xticks( np.arange(30,50,5) )
        plt.axis([xmin,xmax,ymin,ymax])
        plt.xticks( np.arange(xmin,xmax*1.01,5) )
        plt.yticks( np.arange(ymin,ymax*1.01,y_offset) )
        levels = np.linspace(zmin,zmax,51)
        labels = np.linspace(zmin,zmax,6)
        
        plt.title("Time= {0:.0f} min".format(tt/60)) 
            
        plt.contourf(X, Y, Z, levels, cmap=cm.copper, extend="both")
        plt.colorbar(ticks=labels, label='Bed evolution (m)', orientation='horizontal', fraction=0.09, pad=0.3)
        plt.savefig(figname, dpi=600)
        plt.close()
            
    ifile.close()
    
    nn = nn+1

output.close() 
        
        