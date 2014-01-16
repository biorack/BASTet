from omsi.analysis.omsi_analysis_base import omsi_analysis_base
from omsi.analysis.omsi_analysis_data import omsi_analysis_data
from omsi.shared.omsi_dependency import *

class omsi_findpeaks_local(omsi_analysis_base) :
    """Class defining a basic gloabl peak finding. The default implementation computes the peaks on the average
       spectrum and then computes the peak-cube data, i.e., the values for the detected peaks at each pixel.
       
       TODO: The curent version assumes 2D data
    """

    def __init__(self, nameKey="undefined"):
        """Initalize the basic data members"""
        super(omsi_findpeaks_local,self).__init__()
        self.analysis_identifier = nameKey
        self.parameter_names = [ 'msidata' , 'mzdata', 'integration_width', 'peakheight', 'slwindow', 'smoothwidth' , 'printStatus']
        self.data_names = ['peak_mz' , 'peak_value' , 'peak_arrayindex', 'indata_mz']
    
    @classmethod
    def v_qslice(cls , anaObj , z , viewerOption=0) :
        """Implement support for qslice URL requests for the viewer""" 
        #Use the dependency data for slicing here. We do not have a native option to reconstruct images from local peak finding data
        return super(omsi_findpeaks_local,cls).v_qslice( anaObj , z, viewerOption)
    
    @classmethod
    def v_qspectrum( cls, anaObj , x, y , viewerOption=0) :
        """Implement support for qspectrum URL requests for the viewer"""
        #Retrieve the h5py objects for the requried datasets from the local peak finding
        if viewerOption == 0 :
        
            from omsi.shared.omsi_data_selection import check_selection_string, selection_type, selection_to_indexlist
            import numpy as np
            peak_mz = anaObj[ 'peak_mz' ]
            peak_values = anaObj[ 'peak_value' ]
            arrayIndices =  anaObj[ 'peak_arrayindex' ][:]
            indata_mz = anaObj[ 'indata_mz' ]
            #Determine the shape of the original raw data
            if (indata_mz is None) or (arrayIndices is None):
                return None, None
            Nx = arrayIndices[:, 0].max()
            Ny = arrayIndices[:, 1].max()
            NMZ = indata_mz.shape[0]
            numSpectra = arrayIndices.shape[0]
            #Determine the size of the selection and the set of selected items
            xList = selection_to_indexlist(x, Nx)
            yList = selection_to_indexlist(y, Ny)
            if (check_selection_string( x ) == selection_type['indexlist']) and (check_selection_string( y ) == selection_type['indexlist']) :
                if len(xList) == len(yList) :
                    items = [ (xList[i], yList[i]) for i in xrange(0,len(xList))  ]
                else :
                    return None , None
            else :
                items = [0]*(len(xList)*len(yList))
                index = 0
                for xi in xList :
                    for yi in yList :
                        items[index] = (xi , yi )
                        index = index+1
                        
            shapeX = len(items)
            shapeY = 1
            shapeZ = NMZ 
            #Initalize the data cube to be returned
            data = np.zeros( (shapeX, shapeY, shapeZ)  , dtype=peak_values.dtype )
            #Fill the non-zero locations for the data cube with data
            for ni in xrange(0, len(items) ) :
                currentIndex = (items[ni][0]*Ny + items[ni][1])
                currentDX = ni
                currentDY = 0
                startIndex = arrayIndices[ currentIndex ][2]
                if currentIndex < numSpectra :
                    endIndex = arrayIndices[ (currentIndex+1) ][2]
                else :
                    endIndex = peak_values.size
                if startIndex != endIndex :
                    tempValues = peak_values[ startIndex : endIndex ]
                    tempMZ     = peak_mz[ startIndex : endIndex ]
                    data[ currentDX , currentDY , tempMZ ] = tempValues
                else : 
                    #The start and end index may be the same in case that no peaks for found for the given spectrum
                    #The data is already initalized to 0 so there is nothing to do here
                    pass
            
            if len(items) == 1 :
                data = data.reshape( (shapeX , shapeZ) )
              
            #Return the spectra and indicate that no customMZ data values (i.e. None) are needed 
            return data, None
            
        elif viewerOption > 0 :
            return super(omsi_findpeaks_local,cls).v_qspectrum( anaObj , x , y, viewerOption-1)
        else :
            return None, None
        
    @classmethod
    def v_qmz(cls, anaObj, qslice_viewerOption=0, qspectrum_viewerOption=0) :
        """Implement support for qmz URL requests for the viewer"""
        mzSpectra =  None
        labelSpectra = None
        mzSlice = None
        labelSlice = None
        #We do not have native option for qslice, so we rely on the input data in all cases
        if qspectrum_viewerOption == 0 and qslice_viewerOption==0: #Loadings
            mzSpectra = anaObj[ 'indata_mz' ][:]
            labelSpectra = "m/z"
            mzSlice  = None
            labelSlice = None
        elif qspectrum_viewerOption > 0 and qslice_viewerOption>0 :
            mzSpectra, labelSpectra, mzSlice, labelSlice = super(omsi_findpeaks_local,cls).v_qmz( anaObj, qslice_viewerOption=qslice_viewerOption , qspectrum_viewerOption=qspectrum_viewerOption-1)
        elif qspectrum_viewerOption == 0 and qslice_viewerOption>=0 : 
            mzSpectra = anaObj[ 'indata_mz' ][:]
            labelSpectra = "m/z"
            tempA, tempB, mzSlice, labelSlice = super(omsi_findpeaks_local,cls).v_qmz( anaObj, 0 , qspectrum_viewerOption)
            #NOTE: if qspectrum and qslice share the same axis, this call will not return the copied data, i.e., we need to copy the
            #qspectrum values to the qslice values.
            if mzSlice is None :
                mzSlice = tempA
                labelSlice = tempB
        
        return mzSpectra, labelSpectra, mzSlice, labelSlice
        
        
    @classmethod
    def v_qspectrum_viewerOptions(cls , anaObj ) :
        """Define which viewerOptions are supported for qspectrum URL's"""
        dependent_options = super(omsi_findpeaks_local,cls).v_qspectrum_viewerOptions(anaObj)
        re = ["Local peaks"] + dependent_options
        return re

    @classmethod
    def v_qslice_viewerOptions(cls , anaObj ) :
        """Define which viewerOptions are supported for qspectrum URL's"""
        return super(omsi_findpeaks_local,cls).v_qslice_viewerOptions(anaObj)
        
    def execute_analysis(self, msidata , mzdata, integration_width=10, peakheight = 10, slwindow = 100, smoothwidth = 3 , printStatus=False) :
        """Execute the nmf for the given msidata
        
           Keyword Arguments:

           :param msidata: numpy or h5py data object with the msi data. This should be a 3D array. (Mandatory user input)
           :param mzdata: h5py or numpy object with the mzdata for the spectrum. (Mandatory user input)
           :param integration_width: integer paramter indicating the window over which peaks should be integrated. (Default=10)
           :param peakheight: ??? (Default=10)
           :param slwindow: ??? (Default=100)
           :param smoothwidth: ??? (Default=3)
           :param printStatus: Print status messages during execution (Default=False)
           
        """
        #Make sure needed imports are avaiable
        from findpeaks import findpeaks
        import numpy as np
        if printStatus:
            import sys
        
        #Set default parameters if needed
        if not self['integration_width'] :
            self['integration_width'] = 10
        if not self['peakheight'] :
            self['peakheight'] = 10
        if not self['slwindow'] :
            self['slwindow'] = 100
        if not self['smoothwidth'] :
            self['smoothwidth'] = 3
        if not self['printStatus'] :
            self['printStatus'] = False

        #Assign paramters to local variables for convenience
        msidata = self['msidata']
        mzdata = self['mzdata']
        integration_width = self['integration_width']
        peakheight = self['peakheight']
        slwindow = self['slwindow']
        smoothwidth = self['smoothwidth']
        printStatus = self['printStatus']
        
        #Determine the data dimensions
        Nx=msidata.shape[0]
        Ny=msidata.shape[1]
        print msidata.shape
        
        peak_MZ=[]       #The x values for all peaks, stored in a linear array
        peak_values=[]    #The y values for all peaks, stored in a linear array
        peak_arrayindex=np.zeros( shape=(Nx*Ny,3) , dtype='int64' ) #List describing for each pixel the start index where its peaks are stored in the peaks_MZ and peaks_values array
        currentIndex = long(0)
        pixelIndex = 0
        for xi in xrange( 0 , Nx ) :
            for yi in xrange( 0 , Ny ) :
                if printStatus:
                    sys.stdout.write("[" +str( int( 100.* float(pixelIndex)/float(Nx*Ny) )) +"%]"+ "\r")
                    sys.stdout.flush()
                
                
                #Load the spectrum
                y = msidata[xi,yi,:]
                #Find peaks in the spectrum
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
                peak_MZ = peak_MZ + xp
                peak_values = peak_values + yp
                peak_arrayindex[pixelIndex,0] = xi
                peak_arrayindex[pixelIndex,1] = yi
                peak_arrayindex[pixelIndex,2] = currentIndex
                pixelIndex += 1
                currentIndex += len( yp )
        
        #Add the analysis results and parameters to the anlaysis data so that it can be accessed and written to file
        #We here convert the single scalars to 1D numpy arrays to ensure consistency. The data write function can
        #handle also a large range of python built_in types by converting them to numpy for storage in HDF5 but
        #to ensure a consitent behavior we convert the values directly here

        #Save the analysis data to the __data_list so that the data can be saved automatically by the omsi HDF5 file API
        self['peak_mz'] = np.asarray( peak_MZ )
        self['peak_value'] = np.asarray( peak_values )
        self['peak_arrayindex'] = peak_arrayindex
        self['indata_mz'] = mzdata[:]
        

            
def main(argv=None):
    '''Then main function'''

    import sys
    from sys import argv,exit
    from omsi.dataformat.omsi_file import omsi_file
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
    testFPL = omsi_findpeaks_local()
    print "Executing peakfinding analysis"
    testFPL.execute_peakfinding( data , mzdata)
    print "Getting peak finding analysis results"
    pmz = testFPL[ 'peak_mz' ]['data']
    print pmz
    pv = testFPL[ 'peak_value' ]['data']
    print pv
    pai = testFPL[ 'peak_arrayindex' ]['data']
    print pai



if __name__ == "__main__":
    main()

