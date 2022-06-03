#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 26 16:52:11 2022

@authors: Marielle Péré, Asma Chalabi
"""

# ==================================IMPORT===========================================
import numpy as np_
from scipy.linalg import solve
from scipy.linalg import lstsq
import pandas as pd_
import matplotlib.pyplot as plt_
import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
# ==================================IMPORT END===========================================
# ==================================GLOBAL VAR.===========================================
#Marker point coordinates of the FrameSlide in DV
DV_ARRAY=np_.array([[8858.83,	6292.5],
   [-18476.8	,6764.33],
   [-9933.78	,1916.12],
   [9520.25	,-3982.72],
   [-17973.77	,-4115.04],
   [16692.76	,9566.15]])
#Corresponding Marker point coordinates of the FrameSlide in Zeiss
ZEISS_ARRAY=np_.array([[106571.0,	43269.0],
       [106873.0,15901.0],
       [102047.0,24466.0],
       [96294.0	,43917.0],
       [96008.0	,16450.0],
       [109619.0,51127.0]])

#path toward the .pts file of the 20 squares coordinates on the frameSlide for DV
PATH_TO_DV_PTS_FILE=""
#name of .pts file of the 20 squares coordinates on the frameSlide for DV
NAME_DV_PTS_FILE="20170718_points.pts"


#path toward the folder where you want to save .csv file of the 20 squares coordinates converted for Zeiss
PATH_SAVED_ZEISS_FILE=""
#name of the .csv file of the 20 squares coordinates converted for Zeiss
NAME_SAVED_ZEISS_FILE="Zeiss_coordinates_from_DV"
# ===========================FUNCTIONS==================================================

def computeLinearCoefficientDVtoZeiss(
    DV=DV_ARRAY,
    Zeiss=ZEISS_ARRAY):   
    """computeLinearCoefficientDVtoZeiss(
        DV=DV_ARRAY,
        ZEISS=ZEISS_ARRAY)
    Compute the linear coefficients of the linear regression 
    to translate the coordinates used in DeltaVision into the 
    corresponding coordinates on the Zeiss:
        x_Zeiss=a*y_DV+b
        y_Zeiss=c*x_DV+d
        <=> X_Zeiss=A*X_DV+B with X_Zeiss=(x_Zeiss      X_DV=(x_DV
                                           y_Zeiss)           y_DV)
        and A=(0 c
               a 0) 
            B=(b
               d)
    as we rotate the FrameSlide of 90°
    
    Parameters
    -----------------
    DV: numpy array, p x 2, p points on the FrameSlide, 2 coordinates (x,y) in the DeltaVision
        default=DV_ARRAY,
    ZEISS: numpy array, p x 2, p points on the FrameSlide, 2 coordinates (x,y) in the Zeiss
        default=ZEISS_ARRAY
    
    
    
    Returns
    -----------------
    A: numpy array, 2 x 2
    B: numpy array, 2 x 1
    such that X_Zeiss=A x X_DV + B
    
    
    Notes
    -----------------
     
    Example
    -----------------
    see end DVtoZeiss.py
   """
   #init - create the system to solve
    DVcx = np_.ones((DV.shape[0],2))
    DVcx[:,1] = DV[:,0]#[ones(length(DV(:,1)),1) DV(:,1)];
    DVcy = np_.ones((DV.shape[0],2))
    DVcy[:,1] = DV[:,1]
    #solve system
    try:
        regresX = solve(DVcx,Zeiss[:,1])
        regresY = solve(DVcy,Zeiss[:,0])
    except ValueError as ve:
        print("ValueWarning: ", ve," Overdetermined system (several solutions exist) - least square method used instead of solve")
        regresX = lstsq(DVcx,Zeiss[:,1])[0]
        regresY = lstsq(DVcy,Zeiss[:,0])[0]
    # print(regresX)
    # print(regresY)
    #create matrix A and B
    A=np_.array([[0, regresY[1]],[regresX[1],0]])
    B=np_.array([[regresY[0]],[regresX[0]]])
    # % zeiss(:,1)=DV(:,2)*regresY(2)+regresY(1);
    # % zeiss(:,2)=DV(:,1)*regresX(2)+regresX(1);
    return A,B

def computeZeissCoordinatesFromDVcoordinates(
        DV=DV_ARRAY,A=[],B=[],Zeiss_init=[],DV_init=[]):
    """computeZeissCoordinatesFromDVcoordinates(
            DV=DV_ARRAY,A=[],B=[],Zeiss_init=[],DV_init=[])
    
    Compute the linear coefficients of the linear regression 
    to translate the coordinates used in DeltaVision into the 
    corresponding coordinates on the Zeiss by using DV_init and Zeiss_Init
    or, if A and B are defined, use A and B  
    to compute the translation for the DV point in DV:
        x_Zeiss=a*y_DV+b
        y_Zeiss=c*x_DV+d
        <=> X_Zeiss=A*X_DV+B with X_Zeiss=(x_Zeiss      X_DV=(x_DV
                                           y_Zeiss)           y_DV)
        and A=(0 c
               a 0) 
            B=(b
               d)
    as we rotate the FrameSlide of 90°
    
    Parameters
    -----------------
    DV: numpy array, p x 2, p points on the FrameSlide, 2 coordinates (x,y) in the DeltaVision
        default=DV_ARRAY,
    A: list or numpy array, 2 x 2, conversion matrix
        default=[]
        
    B: list or numpy array, 2 x 1, conversion vector
        default=[]
    DV_init: numpy array, m x 2, m marker points on the FrameSlide, 2 coordinates (x,y) in the DeltaVision
        default=[],
    ZEISS_init: numpy array, m x 2, m marker points on the FrameSlide, 2 coordinates (x,y) in the Zeiss
        default=[]
    
    
    
    Returns
    -----------------
    Zeiss_coordinates_array: p x 2, p points, 2 coordinates
    A: numpy array, 2 x 2
    B: numpy array, 2 x 1
    such that X_Zeiss=A x X_DV + B
    
    
    Notes
    -----------------
    If A=[] or B=[], the function compute the translations matrices 
     using Zeiss_init and DV_init, if they are not defined, used the global variable
     
    Example
    -----------------
    see end DVtoZeiss.py
   """
    # if you don't pass the transition matrix A and B, 
    # compute automatically A and B
    # weather with defaults coordinates described in 
    # computeLinearCoefficientDVtoZeiss()
    # or with the coordinates matrices defined in Zeiss_init and DV_init
    
    if len(A)==0 or len(B)==0:
        if Zeiss_init==[]:
            A,B=computeLinearCoefficientDVtoZeiss()
        else:
            A,B=computeLinearCoefficientDVtoZeiss(DV=DV_init,Zeiss=Zeiss_init)
            
    #compute Zeiss coordinated 
    if np_.ndim(DV)==1:
        x=np_.array([[DV[0]],[DV[1]]]) 
        return np_.dot(A,x)+B
    else:
       Zeiss_coordinates_array=np_.zeros(DV.shape)
       for i in range(DV.shape[0]):
           x=np_.array([[DV[i,0]],[DV[i,1]]]) 
           x_Zeiss=np_.dot(A,x)+B
           Zeiss_coordinates_array[i,:]=x_Zeiss.T
       return Zeiss_coordinates_array,A,B


def openDVptsFile(path):
    """takes as input the path to the DV .pts and returns a dataframe
    
    Parameters
    -----------------
    path: str, path to the .pts file
    
    Returns
    -----------------
    df: pandas dataframe, p x 3
        p lines, one line= one point
        3 coordinates: x, y, z
        index: name of the point
        col: [x,y,z]
    
    
    Notes
    -----------------
    
    Example
    -----------------
    see end DVtoZeiss.py"""
    #init
    df=pd_.DataFrame()
    with open(path) as f:
        #read every line of the file
        for rows in f:
            tmp_list=[]
            for chain in rows.split(' '):
                #find the square number
                if ':' in chain:
                    name_col=chain
                #find the x,y,z coordinates
                if '.' in chain :
                    tmp_list.append(float(chain))    
            # print(name_col)
            # print(tmp_list)
            df[name_col]=np_.array(tmp_list)  
    df.index=['x','y','z']
    df=df.T
    return df

def plotFrameSlide(coordinate_array,ax=[]):
    """plotFrameSlide(coordinate_array,ax=[])
    
    plot the Frame slide in 2D with all the point
    in coordinate_array on ax subplot if defined, if not create the figure
    
    Parameters
    -----------------
    coordinate_array: numpy array, p x 2, p point, 2 coordinates
    ax: matplotlib.pyplot AxesSubplot
    
    
    Returns
    -----------------
    ax: matplotlib.pyplot AxesSubplot
    
    Notes
    -----------------
     
    Example
    -----------------
    see end DVtoZeiss.py"""
    
    
    if ax==[]:
        f=plt_.figure()
        ax=f.add_subplot(1,1,1)
    ax.plot(coordinate_array[:,0],coordinate_array[:,1],'ko')
    for i in range(coordinate_array.shape[0]):
        ax.text(coordinate_array[i,0],coordinate_array[i,1],str(i))
    return ax

def plotCoordinatesConversion(DV_array,Zeiss_coordinates_array,title=''):
    """
    plotCoordinatesConversion(DV_array,Zeiss_coordinates_array,title='')
    
    plot the Frame slide in 2D with all the point
    in the DV and in the Zeiss in a 2-subplots figures
    
    Parameters
    -----------------
    DV_array: numpy array, p x 2, p point, 2 coordinates of the points in DV
    Zeiss_array: numpy array, p x 2, p point, 2 coordinates of the points in Zeiss
    title: str, title of the figure
    
    
    Returns
    -----------------
    f: matplotlib.pyplot.figure
    ax1: matplotlib.pyplot AxesSubplot of DV
    ax2: matplotlib.pyplot AxesSubplot of Zeiss
    
    Notes
    -----------------
     
    Example
    -----------------
    see end DVtoZeiss.py"""
    
    f=plt_.figure()
    ax1=f.add_subplot(1,2,1)
    ax1=plotFrameSlide(coordinate_array=DV_array,ax=ax1)
    ax2=f.add_subplot(1,2,2)
    ax2=plotFrameSlide(coordinate_array=Zeiss_coordinates_array,ax=ax2)
    ax1.set_title('DV')
    ax2.set_title('Zeiss')
    f.suptitle(title)
    return f, ax1,ax2

def allConversion(A=[],B=[],Zeiss_init=[],DV_init=[],
                  path_to_dv_pts_file=PATH_TO_DV_PTS_FILE,
                  name_dv_pts_file=NAME_DV_PTS_FILE,
                  path_saved_zeiss_file=PATH_SAVED_ZEISS_FILE,
                  name_saved_zeiss_file=NAME_SAVED_ZEISS_FILE):
    """allConversion(A=[],B=[],path_to_dv_pts_file=PATH_TO_DV_PTS_FILE,
                      name_dv_pts_file=NAME_DV_PTS_FILE,
                      path_saved_zeiss_file=PATH_SAVED_ZEISS_FILE,
                      name_saved_zeiss_file=NAME_SAVED_ZEISS_FILE)
        
    Parameters
    -----------------
    A: list or numpy array, 2 x 2, conversion matrix
        default=[]
        
    B: list or numpy array, 2 x 1, conversion vector
        default=[]
    DV_init: numpy array, m x 2, m marker points on the FrameSlide, 2 coordinates (x,y) in the DeltaVision
        default=[],
    ZEISS_init: numpy array, m x 2, m marker points on the FrameSlide, 2 coordinates (x,y) in the Zeiss
        default=[]
    path_to_dv_pts_file: str, path to DV points .pts file, 
                        default=PATH_TO_DV_PTS_FILE,
    name_dv_pts_file:str, name of the DV .pts file
                    default=NAME_DV_PTS_FILE,
    path_saved_zeiss_file: str,path to the folder where to save the dataframe of the new coordinates
                           default=PATH_SAVED_ZEISS_FILE,
    name_saved_zeiss_file: str, name of the .csv with the new coordinates
                           default=NAME_SAVED_ZEISS_FILE
    
    Returns
    -----------------
    combined_coordinates_system_pts_dataframe: pandas dataframe, p x 4
        p lines, one line= one point
        4 coordinates: x_DV, y_DV, x_Zeiss,y_Zeiss
        index: name of the point
        
    A: numpy array, 2 x 2
    B: numpy array, 2 x 1
    such that X_Zeiss=A x X_DV + B
    
    Notes
    -----------------
    If A=[] or B=[], the function compute the translations matrices 
     using Zeiss_init and DV_init, if they are not defined, used the global variable
     
    Example
    -----------------
    see end DVtoZeiss.py"""
    #load DV points for the 20 squares
    pts_DV_dataframe=openDVptsFile(path_to_dv_pts_file+name_dv_pts_file)
    #turn the dataframe into an array and keep only the x and y coordinates
    DV_array=pts_DV_dataframe[['x','y']].to_numpy()
    #compute Zeiss coordinates from DV coordinates by 
    #automatically compute the conversion matrices
    Zeiss_coordinates_array,A,B=computeZeissCoordinatesFromDVcoordinates(DV=DV_array,\
                                                                        A=A,B=B,Zeiss_init=Zeiss_init,
                                                                        DV_init=DV_init)
    #create the dataframe for Zeiss
    pts_Zeiss_dataframe=pd_.DataFrame(Zeiss_coordinates_array,
                                     index=pts_DV_dataframe.index,
                                     columns=['x_Zeiss','y_Zeiss'])
    
    # create the dataframe on the two coordinate systems together
    pts_DV_dataframe.columns=pts_DV_dataframe.columns+'_DV'
    combined_coordinates_system_pts_dataframe=pd_.concat(( pts_DV_dataframe.iloc[:,:2],pts_Zeiss_dataframe),axis=1)
    
    #save dataframe
    combined_coordinates_system_pts_dataframe.to_csv(path_saved_zeiss_file+name_saved_zeiss_file+".csv",sep='\t')
    
    #plot  
    f,ax1,ax2=plotCoordinatesConversion(DV_array,Zeiss_coordinates_array, title='new test with the  20 squares')
    return combined_coordinates_system_pts_dataframe,A,B
    

if __name__=='__main__': 
    plt_.close('all')
    # ============================================================================= 
    # test: compute translation matrix and then apply the translation 
    # on the DV points obtained with the marker on the FrameSLide
    Zeiss_coordinates_array_test,A_test,B_test=computeZeissCoordinatesFromDVcoordinates()
    # print("True Zeiss coordinates:\n",ZEISS_ARRAY)
    # print("Zeiss coordinates obtained with the translation matrix:\n",Zeiss_coordinates_array_test)
    # f,ax1,ax2=plotCoordinatesConversion(DV_ARRAY,Zeiss_coordinates_array_test,title='Coordinates obtained with the conversion matrices')
    # f,ax1,ax2=plotCoordinatesConversion(DV_ARRAY,ZEISS_ARRAY,title='True coordinates')
    
    # =============================================================================
    #all conversion steps:
    #load DV points for the 20 squares
    pts_DV_dataframe=openDVptsFile(PATH_TO_DV_PTS_FILE+NAME_DV_PTS_FILE)
    #turn the dataframe into an array and keep only the x and y coordinates
    DV_array=pts_DV_dataframe[['x','y']].to_numpy()
    #compute Zeiss coordinates from DV coordinates by 
    #automatically compute the conversion matrices
    Zeiss_coordinates_array,A,B=computeZeissCoordinatesFromDVcoordinates(DV=DV_array)
    #create the dataframe for Zeiss
    pts_Zeiss_dataframe=pd_.DataFrame(Zeiss_coordinates_array,
                                     index=pts_DV_dataframe.index,
                                     columns=['x_Zeiss','y_Zeiss'])
    
    # create the dataframe on the two coordinate systems together
    pts_DV_dataframe.columns=pts_DV_dataframe.columns+'_DV'
    combined_coordinates_system_pts_dataframe=pd_.concat(( pts_DV_dataframe.iloc[:,:2],pts_Zeiss_dataframe),axis=1)
    
    #save dataframe
    combined_coordinates_system_pts_dataframe.to_csv(PATH_SAVED_ZEISS_FILE+NAME_SAVED_ZEISS_FILE+".csv",sep='\t')
    
    #plot  
    f,ax1,ax2=plotCoordinatesConversion(DV_array,Zeiss_coordinates_array, title='new test with the  20 squares')
    # =============================================================================
    #complete pipeline inside one function
    A_final_test=np_.array([[0.        , 10],
           [1.0, 0.        ]])
    B_final_test=np_.array([[1],
           [ 2]])
    combined_coordinates_system_pts_dataframe,A,B=allConversion(A=A_final_test,B=B_final_test)
    combined_coordinates_system_pts_dataframe,A,B=allConversion()
    
    print("coordinates:\n",combined_coordinates_system_pts_dataframe)