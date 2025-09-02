# -*- coding: utf-8 -*-
"""
@Name: Su Hsuan
@Student ID: R12522620
"""

#%% Import libraries
# Additional libraries you might need depending on your task.

# To load an Excel file in Python, you'll need to import the pandas library.
import pandas as pd

# If you're doing numerical operations, you might want to import numpy.
import numpy as np

# If you're working with data visualization, you might want to import matplotlib.
#import matplotlib.pyplot as plt


#%% Load data
df = pd.read_excel('Hw1.xlsx')
data = df.values.tolist()
column_name = df.columns

for i in range(len(column_name)):
    if i == 0 or i%3==1:
        exec(column_name[i]+'=[]')
        for j in range(len(data)):
            if i == 0:
                exec(column_name[i]+'.append(data[j][i])')
            else:
                exec(column_name[i]+'.append(data[j][i:i+3])')


#%% Define Functions
def CoordPelvis(RASI, LASI, RPSI):
    global_coord=[0,0,0]
    X=[1,0,0]
    Y=[0,1,0]
    Z=[0,0,1]
    local_z=np.zeros((len(RASI),3))
    local_y=np.zeros((len(RASI),3))  
    local_x=np.zeros((len(RASI),3))
    Rg2p=np.zeros((len(RASI), 3, 3))
    Vg2p=np.zeros((len(RASI),3))
    RASI_plocal=np.zeros((len(RASI),3))
    LASI_plocal=np.zeros((len(RASI),3))
    RPSI_plocal=np.zeros((len(RASI),3))
    for i in range(len(RASI)):  
        RASI[i] = np.array(RASI[i])
        LASI[i] = np.array(LASI[i])
        RPSI[i] = np.array(RPSI[i])
        local_z[i]=(np.subtract(RASI[i],LASI[i]))/np.linalg.norm(np.subtract(RASI[i],LASI[i]),ord=2)
        local_y[i]=(np.cross(local_z[i],RASI[i]-RPSI[i]))/np.linalg.norm(np.cross(local_z[i],RASI[i]-RPSI[i]),ord=2)
        local_x[i]=np.cross(local_y[i], local_z[i])
        Rg2p[i]=np.array([[np.dot(local_x[i],X),np.dot(local_y[i],X),np.dot(local_z[i],X)],[np.dot(local_x[i],Y),np.dot(local_y[i],Y),np.dot(local_z[i],Y)],[np.dot(local_x[i],Z),np.dot(local_y[i],Z),np.dot(local_z[i],Z)]])
        Vg2p[i]=RASI[i]-global_coord
        RASI_plocal[i]=[0,0,0]
        LASI_plocal[i]= np.dot(np.transpose(Rg2p[i]),(LASI[i]-RASI[i]))
        RPSI_plocal[i]= np.dot(np.transpose(Rg2p[i]),(RPSI[i]-RASI[i]))
       
    return Rg2p, Vg2p, RASI_plocal, LASI_plocal, RPSI_plocal


def CoordThigh(RTRO, RLFC, RMFC):
    global_coord=[0,0,0]
    X=[1,0,0]
    Y=[0,1,0]
    Z=[0,0,1]
    local_z=np.zeros((len(RTRO),3))
    local_y=np.zeros((len(RTRO),3))  
    local_x=np.zeros((len(RTRO),3))
    Rg2t=np.zeros((len(RTRO), 3, 3))
    Vg2t=np.zeros((len(RTRO),3))
    RTRO_tlocal=np.zeros((len(RTRO),3))
    RLFC_tlocal=np.zeros((len(RTRO),3))
    RMFC_tlocal=np.zeros((len(RTRO),3))
    for i in range(len(RTRO)):
        
        
        RTRO[i] = np.array(RTRO[i])
        RLFC[i] = np.array(RLFC[i])
        RMFC[i] = np.array(RMFC[i])
        local_z[i]=(RLFC[i]-RMFC[i])/np.linalg.norm(RLFC[i]-RMFC[i],ord=2)
        local_x[i]=np.cross(RTRO[i]-RLFC[i],local_z[i])/np.linalg.norm(np.cross(RTRO[i]-RLFC[i],local_z[i]),ord=2)
        local_y[i]=np.cross(local_z[i], local_x[i])
        Rg2t[i]=np.array([[np.dot(local_x[i],X),np.dot(local_y[i],X),np.dot(local_z[i],X)],[np.dot(local_x[i],Y),np.dot(local_y[i],Y),np.dot(local_z[i],Y)],[np.dot(local_x[i],Z),np.dot(local_y[i],Z),np.dot(local_z[i],Z)]])
        Vg2t[i]=RTRO[i]-global_coord
        RTRO_tlocal[i]=[0,0,0]
        RLFC_tlocal[i]= np.dot(np.transpose(Rg2t[i]),(RLFC[i]-RTRO[i]))
        RMFC_tlocal[i]= np.dot(np.transpose(Rg2t[i]),(RMFC[i]-RTRO[i]))
        #Rg2t, Vg2t, RTRO_tlocal, RLFC_tlocal, RMFC_tlocal = 0, 0, 0, 0, 0
    
    
    return Rg2t, Vg2t, RTRO_tlocal, RLFC_tlocal, RMFC_tlocal


def CoordShank(RTT, RSHA, RLMA, RMMA):
    global_coord=[0,0,0]
    X=[1,0,0]
    Y=[0,1,0]
    Z=[0,0,1]
    local_z=np.zeros((len(RTT),3))
    local_y=np.zeros((len(RTT),3))  
    local_x=np.zeros((len(RTT),3))
    Rg2s=np.zeros((len(RTT), 3, 3))
    Vg2s=np.zeros((len(RTT),3))
    RTT_slocal=np.zeros((len(RTT),3))
    RSHA_slocal=np.zeros((len(RTT),3))
    RLMA_slocal=np.zeros((len(RTT),3))
    RMMA_slocal=np.zeros((len(RTT),3))
    for i in range(len(RTT)):
        
        RTT[i]=np.array(RTT[i])
        RSHA[i]=np.array(RSHA[i])
        RLMA[i]=np.array(RLMA[i])
        RMMA[i]=np.array(RMMA[i])
        local_x[i]=np.cross(RSHA[i]-RMMA[i],RLMA[i]-RMMA[i])/np.linalg.norm(np.cross(RSHA[i]-RMMA[i],RLMA[i]-RMMA[i]),ord=2)
        local_z[i]=np.cross(local_x[i],(RTT[i]-(RMMA[i]+RLMA[i])/2))/np.linalg.norm(np.cross(local_x[i],(RTT[i]-(RMMA[i]+RLMA[i])/2)),ord=2)
        local_y[i]=np.cross(local_z[i],local_x[i])
        Rg2s[i]=np.array([[np.dot(local_x[i],X),np.dot(local_y[i],X),np.dot(local_z[i],X)],[np.dot(local_x[i],Y),np.dot(local_y[i],Y),np.dot(local_z[i],Y)],[np.dot(local_x[i],Z),np.dot(local_y[i],Z),np.dot(local_z[i],Z)]])
        Vg2s[i]=RTT[i]-global_coord
        RTT_slocal[i]=[0,0,0]
        RSHA_slocal[i]= np.dot(np.transpose(Rg2s[i]),(RSHA[i]-RTT[i]))
        RLMA_slocal[i]= np.dot(np.transpose(Rg2s[i]),(RLMA[i]-RTT[i]))
        RMMA_slocal[i]=np.dot(np.transpose(Rg2s[i]),(RMMA[i]-RTT[i]))
        #Rg2s, Vg2s, RTT_slocal, RSHA_slocal, RLMA_slocal, RMMA_slocal = 0, 0, 0, 0, 0, 0
    
    return Rg2s, Vg2s, RTT_slocal, RSHA_slocal, RLMA_slocal, RMMA_slocal


def CoordFoot(RHEE, RFOO, RTOE):
    global_coord=[0,0,0]
    X=[1,0,0]
    Y=[0,1,0]
    Z=[0,0,1]
    local_z=np.zeros((len(RHEE),3))
    local_y=np.zeros((len(RHEE),3))  
    local_x=np.zeros((len(RHEE),3))
    Rg2f=np.zeros((len(RHEE), 3, 3))
    Vg2f=np.zeros((len(RHEE),3))
    RHEE_flocal=np.zeros((len(RHEE),3))
    RFOO_flocal=np.zeros((len(RHEE),3))
    RTOE_flocal=np.zeros((len(RHEE),3))
    
    for i in range(len(RHEE)):
        
        RHEE[i]=np.array(RHEE[i])
        RFOO[i]=np.array(RFOO[i])
        RTOE[i]=np.array(RTOE[i])
        local_x[i]=(((RFOO[i]+RTOE[i])/2)-RHEE[i])/np.linalg.norm(((RFOO[i]+RTOE[i])/2)-RHEE[i],ord=2)
        local_y[i]=np.cross(local_x[i],RFOO[i]-RTOE[i])/np.linalg.norm(np.cross(local_x[i],RFOO[i]-RTOE[i]),ord=2)
        local_z[i]=np.cross(local_x[i],local_y[i])
        Rg2f[i]=np.array([[np.dot(local_x[i],X),np.dot(local_y[i],X),np.dot(local_z[i],X)],[np.dot(local_x[i],Y),np.dot(local_y[i],Y),np.dot(local_z[i],Y)],[np.dot(local_x[i],Z),np.dot(local_y[i],Z),np.dot(local_z[i],Z)]])
        Vg2f[i]=RHEE[i]-global_coord
        RHEE_flocal[i]=[0,0,0]
        RFOO_flocal[i]= np.dot(np.transpose(Rg2f[i]),(RFOO[i]-RHEE[i]))
        RTOE_flocal[i]= np.dot(np.transpose(Rg2f[i]),(RTOE[i]-RHEE[i]))
        
    #Rg2f, Vg2f, RHEE_flocal, RFOO_flocal, RTOE_flocal = 0, 0, 0, 0, 0
    
    return Rg2f, Vg2f, RHEE_flocal, RFOO_flocal, RTOE_flocal


def CoordG2L(Rg2l, Vg2l, P_global):
    P_local=np.zeros((len(Rg2l),3))   
    for i in range(len(Rg2l)):
        P_local[i]=np.dot(np.transpose(Rg2l[i]),P_global[i]-Vg2l[i])
    return P_local


def CoordL2G(Rg2l, Vg2l, P_local):
    P_global=np.zeros((len(Rg2l),3)) 
    for i in range(len(Rg2l)):
        P_global[i]=np.dot(Rg2l[i],P_local[i])+Vg2l[i]
    return P_global


#Question 1 (40%)
"""
(1) Find the rotation matrix of the pelvis, thigh, shank and
foot coordinate system on the right sides and the corresponding
position vector of the coordinate origin ​​relative to the global
coordinate system over a gait cycle (from frame 237 to frame 450).

(2) Find the local coordinate values of the markers ​​relative to the 
corresponding local coordinate system over a gait cycle (from frame
237 to frame 450).
"""
Rg2p, Vg2p, RASI_plocal, LASI_plocal, RPSI_plocal = CoordPelvis(RASI[236:449], LASI[236:449], RPSI[236:449])
Rg2t, Vg2t, RTRO_tlocal, RLFC_tlocal, RMFC_tlocal = CoordThigh(RTRO[236:449], RLFC[236:449], RMFC[236:449])
Rg2s, Vg2s, RTT_slocal, RSHA_slocal, RLMA_slocal, RMMA_slocal = CoordShank(RTT[236:449], RSHA[236:449], RLMA[236:449], RMMA[236:449])
Rg2f, Vg2f, RHEE_flocal, RFOO_flocal, RTOE_flocal = CoordFoot(RHEE[236:449], RFOO[236:449], RTOE[236:449])
  
    #%% Question 2 (40%)
"""
(1) Determine the motion of the big toe (RBTO) of the right
limb with respect to the foot and shank coordinate system over
a gait cycle (from frame 237 to frame 450).
"""

Rg2l = Rg2s
Vg2l = Vg2s
P_global=RBTO[236:449]
P_local_RBTO_slocal=CoordG2L(Rg2l, Vg2l, P_global)
ans_1=P_local_RBTO_slocal


Rg2l_1 = Rg2f
Vg2l_1 = Vg2f
P_global=RBTO[236:449]
P_local_RBTO_flocal=CoordG2L(Rg2l_1, Vg2l_1, P_global)
ans_2=P_local_RBTO_flocal

"""
Rg2l = Rg2s
Vg2l = Vg2s
P_local=P_local_RBTO_slocal
P_global_from_locals=CoordL2G(Rg2l, Vg2l, P_local)
ans_3=P_global_from_locals

Rg2l_1 = Rg2f
Vg2l_1 = Vg2f
P_local_1=P_local_RBTO_flocal
P_global_from_localf=CoordL2G(Rg2l_1, Vg2l_1, P_local_1)
ans_4=P_global_from_localf
"""
#%% Question 3 & 4 (10%/10%)
"""
Please refer to the attached for more details.
"""