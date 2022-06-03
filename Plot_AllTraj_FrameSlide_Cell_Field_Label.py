"@authors: Asma Chalabi, Marielle Péré "


import numpy as np_
import glob
import pandas as pd_
import matplotlib.pyplot as plt_
from scipy.signal import savgol_filter
import matplotlib

matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42


######## INPUT PARAMETERS ########

# Path to the directory containing .csv files
DIR_PATH = ""

# Time unit to use : "frames", "minutes", "hours"
TIME_UNIT = "minutes"

# Time Between Frames : The acquisition time between two following frames 
TBF = 3 # minutes

# Sliding window for sgolay filter (for trajectories smoothing)
SLD_WIND = 7


def Traj_Proc():
    
    """
    This fuction is used to run the script by calling all the functions.
    
    Input : it doesn't take input parameters.
    
    Output : 
    
    experiment : The experiment organized in a dictionnary.
    time_vec : The time vector.
    labels : Single-cell labels vector.
    mean_derivative : The mean time derivative vector.
    """
    print("\nRead Files ...\n")
    experiment, time_vec, labels = Read_csv_files ()
    
    print("\nProcess Trajectories ...\n")
    experiment = Process_trajectories (experiment)
    
    print("Compute Mean Time Derivative ...\n")
    mean_derivative = Compute_FRET_Derivative (experiment["process. Frame slide"])
    
    print("Sort Cells ...\n")
    sorted_mean_deriv,sorted_labels = SortandSave_FRET_deriv (mean_derivative, labels)
    
    print ("Plot FRET Trajectories ...\n")
    Plot_Top_Bottom_cells (experiment["process. Frame slide"],time_vec, labels,sorted_labels)
    
    print("Done")
    
    
    return experiment, time_vec,labels, mean_derivative




def Read_csv_files () :
    
    """
    This fuction Reads and organizes input files.
    It creates the time vector and the cells' labels.
    
    Input : it doesn't take input parameters.
    
    Output : 
    
    experiment : A dictionnary with separated fields + a matrix with all 
                 merged fields (Frame slide)
                 
    time_vec : The time vector on frame,minutes or hours depending on the 
               user's input.
    
    all_labels : The vector of the cell labels (label = field#.cell#)
    """
    
    path = DIR_PATH
    
    # Add "/" at the end of the path if it doesn't exists
    if path[-1] != "/" :
        path=path+"/"
    
    
    #find all .csv files in the directory
    files = glob.glob(path+"*.csv") 
    
    experiment = {}
    fields_list=[]
    all_labels = []
    
    for file in files :
        
        print(file)
        
        ## Read file
        
        matrix=pd_.read_csv(file, sep="\t")
        matrix=matrix.iloc[:,:].values
    
        trajs_mx =matrix[:,1:] # trajectories matrix
        
        ## Organize the experiment in a dictionary
        
        # Get the field number
        field_num=file.split("/")[-1].split("\\")[1].split(".")[0] # file name : date_siteNum_lineage_condition
        fields_list.append(field_num)
        
        # Add to the dictionnary
        experiment[str(field_num)] = trajs_mx
    
        ## Merge Fields = remodel the Frame Slide
        
        # Add the first field
        if field_num == "01":
            
           experiment["Frame slide"] = experiment[str(field_num)]
        
        # Merge other fields
        else :
            field_1 = experiment["Frame slide"]
            field_2 = trajs_mx
        
            new_shape  = (field_1.shape[0],field_1.shape[1]+field_2.shape[1])
            
            FrameSlide = np_.empty(new_shape)*np_.nan
            FrameSlide[:,:field_1.shape[1]] = field_1
            FrameSlide[:,field_1.shape[1]:new_shape[1]] = field_2 
            experiment["Frame slide"]= FrameSlide
        
        
        ## Create the cell labels
        
        # Call the "Get_cell_labels" function
        labels = Get_cell_labels(field_num,trajs_mx.shape[1])
        all_labels += labels
        
    ## Create time vector
   
    frame_vec = np_.copy(matrix[:,0])
    time_vec = matrix[:,0] 
    time_vec[0] = 0
    
    if TIME_UNIT == "hours":
        time_vec[1:]= (time_vec[1:]-1)*TBF/60
        
    if TIME_UNIT == "minutes":
        time_vec[1:]= (time_vec[1:]-1)*TBF
        
    if TIME_UNIT == "frames":
        time_vec =  frame_vec
        
 
    return experiment, time_vec, np_.asanyarray(all_labels)





def Get_cell_labels (field, cells_num): 
    
    """
    This fuction creates cell labels.
    
    Cell label = field#.cell#
    
    Input :
        
    field : Filed number.
    
    cells_num : Cells number.It conresponds to the column number of field matrix
    
    Output :
    
    labels : A list of labels of all cells in the field
    
    """
    
    labels = []
    
    
    for idx in range(0,cells_num):

        if len(str(idx+1)) > 1  :
            
            label = f"{field}.{str(idx+1)}"
            
        else :
            
            label = f"{field}.0{str(idx+1)}"
        
        labels.append(label)
            
    return labels
            




def Process_trajectories (experiment) :
   
    """
    This fuction processes the FRET trajectories. 
    
    1- Smooth the trajectories with Savitzky-Golay filter.
    2- Normalize the trajectories by subtracting the FRET value at the first frame (t=0) :
       The normalization makes all the trajectories start from 0 in order to compare them.
    
    Input : 
    
    experiment : The experience dictionary.
    
    Output : 
        
    experiment : The experience dictionary with an additionnal matrix of the processed 
                trajectories (process. Frame slide)
    
    """
   
    trajs = experiment["Frame slide"]    # Get all FRET trajectories   
    norm_trajs = np_.empty((trajs.shape))*np_.nan # create an empty matrix
    
    for col in range(0,norm_trajs.shape[1]):
        
        
        traj = trajs[:,col]
    
        # Smooth
        smooth_traj = savgol_filter(traj, SLD_WIND,polyorder=3)
        
       
        # Normalize with the FRET value at the first frame
        norm_traj = smooth_traj- smooth_traj[0]
        
        norm_trajs[:len(smooth_traj),col] = norm_traj
        
    # Save the processed matrix    
    experiment["process. Frame slide"] = norm_trajs
            
    
    return experiment




def Compute_FRET_Derivative (trajs) :
    
    """
    This function computes the mean time derivative of the processed FRET trajectories
    on the last 4th frames.
    
    Input :
        
    trajs : The processed FRET trajectories matrix.(process. Frame slide 
            from the experiment dict).
    
    Output :
    
    mean_derivative : Single-cell mean time derivative vector.
    
    """
 
    fret_derivative = np_.empty((trajs.shape[0],trajs.shape[1]))*np_.NAN # create an empty matrix
    
    
    for traj_idx in range (0, trajs.shape[1]) :
       
        traj = trajs[:,traj_idx]
        
        # Compute the time derivative
        time_derivative = np_.gradient(traj[-5:]) 
        
        fret_derivative[:len(time_derivative),traj_idx] =  time_derivative    
       
    # Compute the mean time derivative
    mean_derivative = np_.nanmean(fret_derivative[:, : ], axis = 0)
                
            
    return mean_derivative




def SortandSave_FRET_deriv (mean_deriv, labels):
    
    """
    This function ranks the mean time derivative values from High to Low values. 
    Then choose the 7 Top High responders and the 7 Bottom Low Responders.
    The result is saved on .xlsx files.
    
    Input :
        
    mean_deriv : Mean time derivative vector.
    
    labels : Cell labels vector.
    
    Output : 
    
    sorted_mean_deriv : Sorted time derivative vector from high to low.
    
    sorted_labels : Sorted labels vector.
    
    dFRET.xlsx : Excel file containing the ranked cells from high to low mean
                time derivative, with their corresponding labels.
                
    dFRET_selected_cells.xlsx : Excel file containing the 7 Top High responders
                ranked from high to low and the 7 Bottom Low Responders ranked 
                from low to high.
    """
    
    # Sort
    sorted_id = np_.argsort(-mean_deriv) #descent
    sorted_mean_deriv = mean_deriv[sorted_id]
    sorted_labels = labels [sorted_id]

    # Save all cells
    out1 = {"dFRET":sorted_mean_deriv, "cell labels": sorted_labels}
    fname1 = "dFRET.xlsx"
    df1 = pd_.DataFrame(out1)
    writer = pd_.ExcelWriter(fname1, engine='xlsxwriter')
    df1.to_excel(writer, index=False, header=["dFRET","Labels"])
    writer.save()
    
    # Save the 7 Top High responders and the 7 Bottom Low Responders
    
    out2 = {"Top dFRET":sorted_mean_deriv[0:7], "Top": sorted_labels[0:7], \
            "Bottom dFRET":sorted_mean_deriv[-7:][::-1],"Bottom":sorted_labels[-7:][::-1]}
    fname2 = "dFRET_selected_cells.xlsx"
    df2 = pd_.DataFrame(out2)
    writer = pd_.ExcelWriter(fname2, engine='xlsxwriter')
    df2.to_excel(writer, index=False, header=["dFRET","Top High Responders","dFRET","Bottom Low Responders"])
    writer.save()
    

    return sorted_mean_deriv,sorted_labels

    

    
def Plot_Top_Bottom_cells (FrameSlide,time, labels,sorted_labels):
    
    """
    This function plots the FRET trajectories.
    The plot highlights the 7 Top High responders (Red) and the 7 Bottom Low 
    Responders (Green) with their corresponding labels.
    
    The labels on the plot correspond to : field#.cell#.rank
    
    Input :
        
    FrameSlide : The processed FRET trajectories matrix.(process. Frame slide 
            from the experiment dict).
    
    time : The time vector.
    
    labels : Cell labels vector.
    
    sorted_labels : Sorted cell labels vector (depending on the mean time derivative
                                               values from high to low)
    
    Output :
        
    FRET.pdf : Plot of the FRET trajectories.
    
    """
    
    # Create figure
    plt_.figure(figsize=(15,8))
    plt_.title("FRET Trajectories")
    plt_.xlim(0, 55)
    
    # Find the top 7 cells
    Top_7 = sorted_labels[0:7]
    top_idx = []
    for top in Top_7 :
        t_idx = np_.where(labels == top)[0][0]
        top_idx.append(t_idx)
    top_traj = FrameSlide[:,top_idx]
    
    # Find the bottom 7 cells
    Bottom_7 = sorted_labels[-7:]
    bottom_idx = []
    for bot in Bottom_7 :
        b_idx = np_.where(labels == bot)[0][0]
        bottom_idx.append(b_idx)
    bottom_traj = FrameSlide[:,bottom_idx]

    # Find the middle cells
    idx_of_interest = top_idx + bottom_idx
    
    # Remove the cell of interest
    middle_cells = np_.delete(FrameSlide,idx_of_interest,1) 
    
    # Plot the middle cells
    plt_.plot(time,middle_cells,"gray")
    

    # Plot the TOP 7 cells in red
    plt_.plot(time,top_traj,"-r",linewidth=3.0)
    
    # Add the TOP 7 cells label
    for top_rank,top_lab in enumerate(Top_7):
       
        plt_.text(time[-1], top_traj[:,top_rank][-1],r"$\leftarrow$"+ f"{top_lab}.{top_rank+1}", fontsize = 6 )
        
    
    # Plot the Bottom 7 cells in green
    plt_.plot(time,bottom_traj, "-g",linewidth=3.0) 
    
    # Add the Bottom 7 cells label
    incr_rank = 7
    for bot_rank,bot_lab in enumerate(Bottom_7):
       
        plt_.text(time[-1], bottom_traj[:,bot_rank][-1],r"$\leftarrow$"+ f"{bot_lab}.{incr_rank}", fontsize = 6)
        incr_rank-=1
    
    # Save figure
    plt_.savefig("FRET.pdf",transparent= True)



###### MAIN ###### 


experiment, time, labels, mean_deriv = Traj_Proc()









