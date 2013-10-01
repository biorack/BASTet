from omsi.analysis.omsi_analysis_base import omsi_analysis_base
from omsi.analysis.omsi_analysis_data import omsi_analysis_data
from omsi.shared.omsi_dependency import *

class omsi_findpeaks_global(omsi_analysis_base) :
    """Class defining a basic gloabl peak finding. The default implementation computes the peaks on the average
       spectrum and then computes the peak-cube data, i.e., the values for the detected peaks at each pixel.
       
       TODO: The curent version assumes 2D data
    """

    def __init__(self, nameKey="undefined"):
        """Initalize the basic data members"""
        super(omsi_findpeaks_global,self).__init__()
        self.analysis_identifier = nameKey
        
    @classmethod
    def v_qslice(cls , anaObj , z , viewerOption=0) :
        """Implement support for qslice URL requests for the viewer""" 
        if viewerOption == 0 :
            dataset =  anaObj[ 'peak_cube' ]
            try:
                data = eval("dataset[:,:, %s]" %(z,))
                return data
            except:
                print "Global peak selection failed"
                return None
        elif viewerOption >= 0 :
            return super(omsi_findpeaks_global,cls).v_qslice( anaObj , z, viewerOption-1)
        else :
            return None
    @classmethod
    def v_qspectrum( cls, anaObj , x, y , viewerOption=0) :
        """Implement support for qspectrum URL requests for the viewer"""
        #Get the h5py dataset with the peak_cube data
        data = None
        customMZ = None
        if viewerOption == 0 :
            from omsi.shared.omsi_data_selection import check_selection_string, selection_type
            dataset =  anaObj[ 'peak_cube' ]
            if (check_selection_string(x) == selection_type['indexlist']) and (check_selection_string(y) == selection_type['indexlist']) :
                #The peak-cube data is usually small enough. To handle the multiple list selection case
                #we here just load the full data cube and use numpy to do the subselection. Note, this
                #version would work for all selection types but we would like to avoid loading the 
                #full data if we don't have to.
                data = eval("dataset[:][%s,%s, :]" %(x,y))
            else :
                data = eval("dataset[%s,%s, :]" %(x,y))
            #Return the spectra and indicate that no customMZ data values (i.e. None) are needed 
            return data, None
        elif viewerOption > 0 :
            return super(omsi_findpeaks_global,cls).v_qspectrum( anaObj , x , y, viewerOption-1)
            
        return data , customMZ
        
    @classmethod
    def v_qmz(cls, anaObj, qslice_viewerOption=0, qspectrum_viewerOption=0) :
        """Implement support for qmz URL requests for the viewer"""
        mzSpectra =  None
        labelSpectra = None
        mzSlice = None
        labelSlice = None
        #We do not need to handle the qslice_viewerOption separately here since there is only one option right now
        if qspectrum_viewerOption == 0 and qslice_viewerOption==0: #Loadings
            mzSpectra =  anaObj[ 'peak_mz' ][:]
            labelSpectra = "m/z"
            mzSlice  = None
            labelSlice = None
        elif qspectrum_viewerOption > 0 and qslice_viewerOption>0 :
            mzSpectra, labelSpectra, mzSlice, labelSlice = super(omsi_findpeaks_global,cls).v_qmz( anaObj, qslice_viewerOption-1 , qspectrum_viewerOption-1)
        elif qspectrum_viewerOption == 0 and qslice_viewerOption>0 :
            mzSpectra =  anaObj[ 'peak_mz' ][:]
            labelSpectra = "m/z"
            tempA, tempB, mzSlice, labelSlice = super(omsi_findpeaks_global,cls).v_qmz( anaObj, qslice_viewerOption-1 , 0)
        elif qspectrum_viewerOption > 0 and qslice_viewerOption==0 :
            mzSlice =  anaObj[ 'peak_mz' ][:]
            labelSlice = "m/z"
            mzSpectra, labelSpectra, tempA, tempB = super(omsi_findpeaks_global,cls).v_qmz( anaObj, 0 , qspectrum_viewerOption-1)
        
        return mzSpectra, labelSpectra, mzSlice, labelSlice
        
    @classmethod
    def v_qspectrum_viewerOptions(cls , anaObj ) :
        """Define which viewerOptions are supported for qspectrum URL's"""
        dependent_options = super(omsi_findpeaks_global,cls).v_qspectrum_viewerOptions(anaObj)
        re = ["Peak cube"] + dependent_options
        return re

    @classmethod
    def v_qslice_viewerOptions(cls , anaObj ) :
        """Define which viewerOptions are supported for qspectrum URL's"""
        dependent_options = super(omsi_findpeaks_global,cls).v_qslice_viewerOptions(anaObj)
        re = ["Peak cube"] + dependent_options
        return re

    
    def execute_peakfinding(self, msidata , mzdata, integration_width=10, peakheight = 2, slwindow = 100, smoothwidth = 3, msidata_dependency=None) :
        """Execute the global peak finding for the given msidata
        
           Keyword Arguments:

           :param msidata: numpy or h5py data object with the msi data. This should be a 3D array.
           :param mzdata: h5py or numpy object with the mzdata for the spectrum
           :param integration_width: integer paramter indicating the window over which peaks should be integrated
           :param peakheight: ???
           :param slwindow: ???
           :param smoothwidth: ???
           :param msidata_dependency: The path or omsi_file_ object on which the msidata input data object depends/originates from. 

        """
        from findpeaks import findpeaks
        import numpy as np
        
        #Load the input data
        data = msidata[:]
        #Determine the data dimensions
        Nx=data.shape[0]
        Ny=data.shape[1]
        Nz=data.shape[2]
        #Linearize the spectral data
        y = data.reshape(Nx*Ny,Nz)
        #Compute the average spectrum
        y = np.mean(y,axis=0)
        #Find peaks in the average spectrum
        A = findpeaks(mzdata[:], y, smoothwidth, slwindow, peakheight)
        y = A.smoothListGaussian()
    
        # from the smoothed spectra subtract a sliding minima 
        A = findpeaks(mzdata[:],y, smoothwidth, slwindow, peakheight)
        slmin = [ x for x in A.sliding_window_minimum()]
        y = y - slmin
        
        # find peaks in the smoothed, background subtracted spectra
        A = findpeaks(mzdata[:],y, smoothwidth, slwindow, peakheight)
        [pkmax,pkmin]=A.peakdet()
        xp = [ x[0] for x in pkmax  ]
        yp = [ x[1] for x in pkmax  ]
        
        # integrate peaks +/- integration_width bins around each of the peaks found in the total spectra
        # TODO : THIS LOOP NEEDS TO BE CONVERTED TO A MAX INSTEAD OF A SUM
        im = data[:,:,xp]
        im = im*0.0;
        xp = np.array(xp)
        for i in range(-integration_width,integration_width):
            idx = xp+i
            im = im + data[:,:,idx]
        
        #Add the analysis results and parameters to the anlaysis data so that it can be accessed and written to file
        #We here convert the single scalars to 1D numpy arrays to ensure consistency. The data write function can
        #handle also a large range of python built_in types by converting them to numpy for storage in HDF5 but
        #to ensure a consitent behavior we convert the values directly here
            
        #Clear any previously stored analysis data
        self.clear_analysis_data()
        self.clear_parameter_data()
        self.clear_dependency_data()
        
        #Save the analysis data to the __data_list so that the data can be saved automatically by the omsi HDF5 file API
        self.add_analysis_data( name='peak_cube' , data=im , dtype=str(im.dtype) ) 
        self.add_analysis_data( name='peak_mz' , data=mzdata[xp] , dtype=str(mzdata.dtype) ) 
        
        #Save the analysis parameters to the __parameter_list so that the data can be saved automatically by the omsi HDF5 file API
        iw = np.asarray( [ integration_width ] )
        self.add_parameter_data( name='integration_width' , data=iw , dtype=str(iw.dtype) ) 
        ph = np.asarray( [ peakheight ] )
        self.add_parameter_data( name='peakheight' , data=ph , dtype=str(ph.dtype) ) 
        slw = np.asarray( [ slwindow ] )
        self.add_parameter_data( name='slwindow' , data=slw , dtype=str(slw.dtype) ) 
        smw = np.asarray( [ smoothwidth ] )
        self.add_parameter_data( name='smoothwidth ' , data=smw  , dtype=str(smw.dtype) ) 
        
        #Save the analysis dependencies to the __dependency_list so that the data can be saved automatically by the omsi HDF5 file API
        if msidata_dependency is not None :
            if isinstance( msidata_dependency , omsi_dependency ) :
                #Add the depency as given by the user
                self.add_dependency_data( msidata_dependency )
            else :
                #The user only gave us the object that we depend on so we need to construct the 
                self.add_dependency_data( omsi_dependency( param_name = 'msidata', link_name='msidata', omsi_object=msidata_dependency, selection=None ) )


def main(argv=None):
    '''Then main function'''

    import sys
    from sys import argv,exit
    from omsi.dataformat.omsi_file import omsi_file
    import matplotlib.pyplot as plt
    import matplotlib.gridspec as gridspec
    import numpy as np

    if argv is None:
        argv = sys.argv

    #Check for correct usage
    if len(argv) <2 :
        print "USAGE: Call \"omsi_findpeaks_global OMSI_FILE [expIndex dataIndex]   \" "
        print "\n"
        print "This is a simple test function to test the peak finding on a given omsi HDF5 file"
        exit(0)

    #Read the input arguments 
    omsiInFile  = argv[1]
    expIndex = 0
    dataIndex = 0
    if len(argv)==4 :
        expIndex = int(argv[2] )
        dataIndex = int(argv[3] )

    #Open the input HDF5 file
    try:
        omsiFile = omsi_file( omsiInFile , 'r' ) #Open file in read only mode
    except:
        print "Unexpected openeing the input  file:", sys.exc_info()[0]
        exit(0)
        
    #Get the experiment and dataset
    exp = omsiFile.get_exp( expIndex )
    data = exp.get_msidata( dataIndex )
    mzdata = exp.get_instrument_info().get_instrument_mz()

    #Execute the peak finding
    testFPG = omsi_findpeaks_global()
    print "Executing peakfinding analysis"
    testFPG.execute_peakfinding( data , mzdata)
    print "Getting peak finding analysis results"
    peakCube = testFPG[ 'peak_cube' ]['data']
    print peakCube
    
    #Plot the first three peak images
    print "Plotting example peak finding analysis results"
    Nx = peakCube.shape[0]
    Ny = peakCube.shape[1]
    ho1 = peakCube[:,:,0];
    ho2 = peakCube[:,:,0];
    ho3 = peakCube[:,:,0];
    
    mainFig = plt.figure()
    gs = gridspec.GridSpec(1, 4)
    imageFig = mainFig.add_subplot(gs[0])
    imageFig.autoscale(True,'both',tight=True)
    imagePlot = imageFig.pcolor( np.log (ho1 + 1) ) 

    imageFig = mainFig.add_subplot(gs[1])
    imageFig.autoscale(True,'both',tight=True)
    imagePlot = imageFig.pcolor( np.log (ho2 + 1) )

    imageFig = mainFig.add_subplot(gs[2])
    imageFig.autoscale(True,'both',tight=True)
    imagePlot = imageFig.pcolor( np.log (ho3 + 1) )


    #do the three color
    ho = np.zeros( shape=(Nx, Ny , 3 ) )
    temp = np.log ( peakCube[:,:,0] + 1)
    temp = temp - temp.min()
    temp = temp / temp.max()
    ho[:,:,0] = temp
    
    temp = np.log ( peakCube[:,:,1] + 1)
    temp = temp - temp.min()
    temp = temp / temp.max()
    ho[:,:,1] = temp
    
    temp = np.log ( peakCube[:,:,2] + 1)
    temp = temp - temp.min()
    temp = temp / temp.max()
    ho[:,:,2] = temp
    
    imageFig = mainFig.add_subplot(gs[3])
    imageFig.autoscale(True,'both',tight=True)
    imagePlot = imageFig.imshow( ho )
   
    plt.show()
    
    

   
if __name__ == "__main__":
    main()

