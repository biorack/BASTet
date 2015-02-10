from omsi.dataformat.omsi_file import *
#Open the file
#f is an object of the tyep omsi_file. For more information execute:
#help( omsi_file )
f = omsi_file("/project/projectdirs/openmsi/omsi_data/20120711_Brain.h5" , 'r' )
#Get the number of experiments
num_exp = f.get_num_experiments()
#Get the first experiment. 
#exp0  is an object of the type omsi_file_experiment For more information execute:
#help( omsi_file_experiment )
exp0 = f.get_experiment(0)
#Get the number of available datasets and analyses
num_dat = exp0.get_num_msidata()
num_ana = exp0.get_num_analysis()
#Get the raw MSI dataset, i.e., the API object that manages the dataset
#We can slice into d using standard array syntax. The MSI dataset defines
#a 3D cube in (x,y,m/z) and contains the corresponding intensity values.
#The m/z values are the same for each spectrum and are stored in a 
#separate 1D array. In the call below we get the experiment with index 0
#and the msi dataset with index 0 associated with the experiment. An
#OpenMSI file may contain multipe experiments, each of which may have
#multiple MSI datasets associated with it. For now, all data you have.
#d is an object of the type omsi_file_msidata. For more information execute:
#help( omsi_file_msidata )
d = exp0.get_msidata(0)
#Number of pixel in x, y , m/z
numX  = d.shape[0]
numY  = d.shape[1]
numMZ = d.shape[2]
#Get the m/z data. This is a 1D array of length d.shape[2]
mz = d.mz
#Load the 10th ion image
image10 = d[:,:,10]
#Load spectrum of pixel (2,3) 
spectrum_2_3 = d[2,3,:]



