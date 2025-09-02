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
import matplotlib.pyplot as plt
import sympy 
import math

#%% Load data
df = pd.read_excel('Hw3.xlsx')
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

# Hw1
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

# Hw2

def RotFormula(sequence):
    t1, t2, t3 = sympy.symbols('t1 t2 t3')
    R = np.eye(3)
    for axis, angle_symbol in zip(sequence, [t1, t2, t3]):
        if axis == 'x':
            rotation_matrix = np.array([[1, 0, 0],[0, sympy.cos(angle_symbol), -sympy.sin(angle_symbol)],[0, sympy.sin(angle_symbol), sympy.cos(angle_symbol)]])
        elif axis == 'y':
            rotation_matrix = np.array([[sympy.cos(angle_symbol), 0, sympy.sin(angle_symbol)],[0, 1, 0],[-sympy.sin(angle_symbol), 0, sympy.cos(angle_symbol)]])
        elif axis == 'z':
            rotation_matrix = np.array([[sympy.cos(angle_symbol), -sympy.sin(angle_symbol), 0],[sympy.sin(angle_symbol), sympy.cos(angle_symbol), 0],[0, 0, 1]])

        R = np.dot(R, rotation_matrix)

    return R



def Rot2Ang(Rot, sequence):
    
    theta = np.zeros((len(Rot),3))
    t1=np.zeros(len(Rot))
    t2=np.zeros(len(Rot))
    t3=np.zeros(len(Rot))
    for i in range (len(Rot)):
        R_3x3 = Rot[i]       
        t2[i]=np.arcsin(R_3x3[2,1])
        t1[i]=np.arctan2(-R_3x3[0,1],R_3x3[1,1])
        t3[i]=np.arctan2(-R_3x3[2,0],R_3x3[2,2])
        
        theta[i] = [t1[i], t2[i] ,t3[i]]            

    return theta


def Ang2Rot(theta, sequence):
    Rot = np.eye(3)
    for axis, angle in zip(sequence, theta):
        if axis == 'x':
            rotation_matrix = np.array([[1, 0, 0],
                                        [0, np.cos(angle), -np.sin(angle)],
                                        [0, np.sin(angle), np.cos(angle)]])
        elif axis == 'y':
            rotation_matrix = np.array([[np.cos(angle), 0, np.sin(angle)],
                                        [0, 1, 0],
                                        [-np.sin(angle), 0, np.cos(angle)]])
        elif axis == 'z':
            rotation_matrix = np.array([[np.cos(angle), -np.sin(angle), 0],
                                        [np.sin(angle), np.cos(angle), 0],
                                        [0, 0, 1]])

        Rot = np.dot(Rot, rotation_matrix)
    return Rot 


def RotAngConvert(input, sequence):
    if len(input.shape) == 3:
        
        return Rot2Ang(input, sequence)
    elif len(input.shape) == 1:
       
        return Ang2Rot(input, sequence)



#%% Hw3
def PolyFit(xi, yi, n):
    Y=np.zeros((len(yi),1))
    X=np.zeros((len(xi),n+1))
    p=np.zeros((n+1,1))
    for row in range(len(xi)):
        for col in range(n+1):
            X[row, col] = xi[row] ** (n - col)
            Y[row] = yi[row]  
    p = np.linalg.inv(np.transpose(X) @ X) @ np.transpose(X) @ Y
    
    return p


def PolyDer(p, dorder):
    n = len(p) - 1
    dp = np.zeros((n - dorder + 1, 1))
    for i in range(n - dorder + 1):
        dp[i, 0] = (np.math.factorial(n - i) / np.math.factorial(n - i - dorder)) * p[i, 0]
    return dp


def PolyVal(p, xi):
    n = len(p) - 1
    yi = np.zeros(len(xi))
    for i in range(len(xi)):
        for j in range(n + 1):
            yi[i] += p[j] * (xi[i] ** (n - j))
    return yi



def Derivative(xi, yi , dorder):
    p = PolyFit(xi, yi, 5)
    dp = PolyDer(p, dorder)
    dyi = PolyVal(dp, xi)
    return dyi



def Ang2LocalAngular(theta, seq, smprate):
    nframes = len(theta)
    AngVel = np.zeros((nframes, 3))
    AngAcc = np.zeros((nframes, 3))

    for i in range(nframes):
        R = RotFormula(seq)(*theta[i])
        omega = np.zeros((3, 3))
        alpha = np.zeros((3, 3))

        for j in range(3):
            axis = seq[j]
            axis_index = ord(axis) - ord('x')
            omega[axis_index, j] = 1

        AngVel[i] = np.dot(omega, theta[i])
        alpha = Derivative([smprate * k for k in range(nframes)], AngVel[:, j], 1)

        AngAcc[i] = np.dot(alpha, theta[i])

    return AngVel, AngAcc
   


def Rot2LocalAngular(Rg2l, smprate):
    nframes = len(Rg2l)
    AngVel = np.zeros((nframes, 3))
    AngAcc = np.zeros((nframes, 3))

    for i in range(nframes):
        R = Rg2l[i]
        omega = np.zeros((3, 3))
        alpha = np.zeros((3, 3))

        for j in range(3):
            omega[j, :] = Derivative([smprate * k for k in range(nframes)], R[:, j], 1)
            alpha[j, :] = Derivative([smprate * k for k in range(nframes)], omega[j, :], 1)

        AngVel[i] = np.dot(np.linalg.inv(R), omega)
        AngAcc[i] = np.dot(np.linalg.inv(R), alpha)

    return AngVel, AngAcc


def Hw3_Q1_f(x):
    # Define the function for Hw3's Q1
    return 1/(8*x**2-16*x+9)


def Hw3_Q1_df(x):
    # Define the function for Hw3's Q1
    return -(16*x-16)/(8*x**2-16*x+9)**2


def Hw3_Q1_visualisation(degrees: list, point_numbers: list):
    if not(len(degrees)==len(point_numbers)):
        print("Error:")
        print("The number of elements in the 'degrees' list must match the number of elements in the 'point_numbers' list.")
        return
    
    # Generate equally spaced points
    x_interp = np.linspace(0, 2, 10000)
    plt.figure(figsize=(12, 8))
    plt.rcParams['figure.dpi'] = 300
    fig, axs = plt.subplots(2, 2, figsize=(10, 4))
    
    # Interpolate using polynomials of different degrees
    for i in range(len(degrees)):
        degree = degrees[i]
        point_number = point_numbers[i]
        x = np.linspace(0, 2, point_number)
        y = Hw3_Q1_f(x)

        # Perform polynomial interpolation
        x = x.tolist()
        y = y.tolist()
        
        coefficients_0dorder = np.array(PolyFit(x, y, degree))  # TODO
        coefficients_1dorder = np.array(PolyDer(coefficients_0dorder, 1))  # TODO

        y_interp_0dorder = np.array(PolyVal(coefficients_0dorder, x_interp))  # TODO
        y_interp_1dorder = np.array(PolyVal(coefficients_1dorder, x_interp))  # TODO
        
        
        # Plot the interpolation
        axs[1,0].plot(x_interp, y_interp_0dorder, label=f'Degree: {degree}; Point Number: {point_number}')
        axs[1,1].plot(x_interp, y_interp_1dorder, label=f'Degree: {degree}; Point Number: {point_number}')

    # Plot the original function for reference
    y_true = Hw3_Q1_f(x_interp)
    dy_true = Hw3_Q1_df(x_interp)
    
    for i in range(2):
        axs[i,0].plot(x_interp, y_true, label='True Function', linewidth=2, linestyle=':', color='k')
        axs[i,1].plot(x_interp, dy_true, label='True Function', linewidth=2, linestyle=':', color='k')
        for j in range(2):
            axs[i,j].grid(True)

    axs[0,0].set_title("Results for f(x) = 1/(8*x^2-16*x+9)")
    axs[0,1].set_title("Results for f'(x) = -(16*x-16)/(8*x^2-16*x+9)^2")
    axs[1,1].legend(bbox_to_anchor=(1.05, 1.0), loc='upper left')
    
    plt.show()


#%% Question 1 (30%)
"""
(1) Compare the values of the original function with those obtained
 from the fitted polynomial equation to investigate how different
 polynomial orders and sampling point numbers affect the fitting
 results of the function.
"""
Hw3_Q1_visualisation([7,10,13], [15,15,15])
Hw3_Q1_visualisation([13,13,13], [15,30,45])


#%% Question 2 (10%)

"""
(1) Please utilise the least squares method to determine the
 coefficients of the fifth-order polynomial equation that fits
 the curve of knee joint angle in the sagittal plane over a gait
 cycle.

(2) Please calculate the angular velocity and angular acceleration
 of the knee joint in the sagittal plane.
"""
Rg2t, Vg2t, RTRO_tlocal, RLFC_tlocal, RMFC_tlocal = CoordThigh(RTRO[236:449], RLFC[236:449], RMFC[236:449])
Rg2s, Vg2s, RTT_slocal, RSHA_slocal, RLMA_slocal, RMMA_slocal = CoordShank(RTT[236:449], RSHA[236:449], RLMA[236:449], RMMA[236:449])
rRt2s = np.zeros((len(Rg2t),3,3))


for i in range (len(Rg2t)):
    rRt2s[i] = np.dot(np.transpose(Rg2t[i]),Rg2s[i])
    sequence = 'zxy'
    result2 = Rot2Ang(rRt2s,sequence)
Time_array = np.array(Time[236:449])
yi_Q2 = result2[:,0]
#yi_Q2 = np.degrees(result2[:,0])
xi_Q2 = Time_array
n = 5
p_Q2 = PolyFit(xi_Q2,yi_Q2,n)
ANS2_1 = Derivative(xi_Q2, yi_Q2 , 1)
ANS2_2 = Derivative(xi_Q2, yi_Q2 , 2)


#%% Question 3 (40%)
"""
(1) Please calculate the pelvic angular velocity and angular
 acceleration relative to its own local coordinate system during
 a gait cycle. 
(2) Compare the pelvic angular acceleration
 (i) obtained by using the reference formula and
 (ii) obtained by numerically differentiating the angular velocity from discrete data values.
 #using function:
 Rg2p, Vg2p, RASI_plocal, LASI_plocal, RPSI_plocal = CoordPelvis(RASI[236:449], LASI[236:449], RPSI[236:449])
 theta = np.degrees(Rot2Ang(Rg2p,'zxy'))
 seq='zxy'
 smprate=np.array(Time)
 Ans3_omega,ANS3_alpha = Ang2LocalAngular(theta, seq, smprate)
"""

#計算角速度與角加速度
Rg2p, Vg2p, RASI_plocal, LASI_plocal, RPSI_plocal = CoordPelvis(RASI[236:449], LASI[236:449], RPSI[236:449])
#theta = np.degrees(Rot2Ang(Rg2p,'zxy'))
theta = Rot2Ang(Rg2p,'zxy')
z_ang=theta[:,0]
x_ang=theta[:,1]
y_ang=theta[:,2]
der_z = Derivative(Time_array,theta[:,0],1)
der_x = Derivative(Time_array,theta[:,1],1)
der_y = Derivative(Time_array,theta[:,2],1)

der2_z = Derivative(Time_array, der_z, 1)
der2_x = Derivative(Time_array, der_x, 1)
der2_y = Derivative(Time_array, der_y, 1)

ANS3_omega = np.zeros((len(theta),3))
ANS3_alpha = np.zeros((len(theta),3))
for i in range(len(theta)):
    X_rot = np.array([[1, 0, 0],
                            [0, np.cos(x_ang[i]), -np.sin(x_ang[i])],
                            [0, np.sin(x_ang[i]), np.cos(x_ang[i])]])
    Y_rot = np.array([[np.cos(y_ang[i]), 0, np.sin(y_ang[i])],
                             [0, 1, 0],
                             [-np.sin(y_ang[i]), 0, np.cos(y_ang[i])]])
    Z_rot = np.array([[np.cos(z_ang[i]), -np.sin(z_ang[i]), 0],
                            [np.sin(z_ang[i]), np.cos(z_ang[i]), 0],
                            [0, 0, 1]])

    ANS3_omega[i] = Y_rot @ X_rot @ np.array([0, 0, der_z[i]]) + Y_rot @ np.array([der_x[i], 0, 0]) + np.array([0, der_y[i], 0])   
    ANS3_alpha[i] = Y_rot @ X_rot @ np.array([0, 0, der2_z[i]]) + Y_rot @ np.array([der2_x[i], 0, 0]) + np.array([0, der2_y[i], 0])  


#以參考公式計算骨盆的角加速度
Time_array_Q3 = np.array(Time[236:449])
alpha_ref = np.zeros((len(theta), 3))
for i in range(len(theta)):
    alpha_ref[i] = np.array([
        [-der_y[i]*np.sin(y_ang[i]), 0, -der_x[i]*np.sin(x_ang[i])*np.sin(y_ang[i])+der_y[i]*np.cos(x_ang[i])*np.cos(y_ang[i])],
        [0, 0, -der_x[i]*np.cos(x_ang[i])],
        [-der_y[i]*np.cos(y_ang[i]), 0, -der_x[i]*np.sin(x_ang[i])*np.cos(y_ang[i])-der_y[i]*np.cos(x_ang[i])*np.sin(y_ang[i])]
    ]) @ np.array([der_x[i], der_y[i], der_z[i]]) + np.array([
        [np.cos(y_ang[i]), 0, np.cos(x_ang[i]*np.sin(y_ang[i]))],
        [0, 1, -np.sin(x_ang[i])],
        [-np.sin(y_ang[i]), 0, np.cos(x_ang[i])*np.cos(y_ang[i])]
    ]) @ np.array([der2_x[i], der2_y[i], der2_z[i]])
        
#繪圖比較
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 10))
ax1.plot(Time_array_Q3, ANS3_alpha[:, 0], label='Z-axis(calculated)')
ax1.plot(Time_array_Q3, ANS3_alpha[:, 1], label='X-axis(calculated)')
ax1.plot(Time_array_Q3, ANS3_alpha[:, 2], label='Y-axis(calculated)')
ax1.legend()
ax1.set_title('Angular Acceleration Over Time(calculated)')
ax1.set_xlabel('Time')
ax1.set_ylabel('Angular Acceleration')
ax1.grid(True)

ax2.plot(Time_array, alpha_ref[:, 0], label='Z-axis (Reference Formula)')
ax2.plot(Time_array, alpha_ref[:, 1], label='X-axis (Reference Formula)')
ax2.plot(Time_array, alpha_ref[:, 2], label='Y-axis (Reference Formula)')
ax2.legend()
ax2.set_title('Angular Acceleration Over Time(reference formula)')
ax2.set_xlabel('Time')
ax2.set_ylabel('Angular Acceleration')
ax2.grid(True)
plt.show()


#%% Question 3 & 4 (10%/10%)
"""
Please refer to the attached for more details.
"""