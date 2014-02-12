import numpy as np
import sys, os, time


#from omsi.dataformat.omsi_file import omsi_file, omsi_file_analysis
#from omsi.analysis.omsi_analysis_base import omsi_analysis_base

"""This an extended list of types indicated by the check_selection_string function. Indicies <0 are assumed to be invalid selections."""
selection_type = {'invalid':-1, 'index':0, 'indexlist':2, 'all':3, 'range':4}


def check_selection_string(selectionString):
    """Check whether the given selection string is valid, and indicate which type of selection
       the string defined. Checking the selection string is meant as a safeguard to prevent 
       attakers from being able to insert malitious code.
       
       :param selectionString: String given by the user with the desired selection
       
       :returns: String indicating the type of selection as defined in selection_type:
       
            * 'indexlist' : Selection of the form [1,2,3]
            * 'all' : Selection of the form ':'
            * 'range: Selection of the form 'a:b'
            * 'index: A single index selection, e.g., '1'
            * 'invalid': An unsupported selection 
                 
       """
    import re
    #Check if we actually have a string
    if selectionString is None:
        return selection_type['invalid']
    #Check if the selection defined a list of indicies
    elif re.match('^\[([0-9]+,)*[0-9]+\]$', selectionString): 
        return selection_type['indexlist'] 
    #Check if the selection defines and ":" all selection
    elif re.match('^:$', selectionString): 
        return selection_type['all']
    #Check if the selection defines a range  a:b type selection
    elif re.match('^[0-9]+:-?[0-9]+$', selectionString):   
        return selection_type['range']
    #Check if the selection defines a single index value type selection 
    elif re.match('^-?[0-9]+$', selectionString):
        return selection_type['index']
    #If none of the above selection types are given, then the selection is declared invalid 
    else:
        return selection_type['invalid']

def selection_string_to_object(selectionString):
    """Convert the given selection string to a python selection object, i.e., either a slice, list or integer index.
    
       :param selectionString: A selection string of the type indexlist
       
       :returns: 
       
            * An integer index if an index selection is specified
            * A python list of indices if a list specified in the string
            * A python slice opject if a slice opteration is specified by the string
            
    """
    selectionType = check_selection_string(selectionString) 
    if selectionType == selection_type['indexlist']:
        stringList = selectionString[1:-1].split(",")
        try:
            parsedList = [int(i) for i in stringList]
        except:
            return None
        return parsedList
    elif selectionType == selection_type['index']:
        try: 
            return int(selectionString)
        except:
            return None
    elif selectionType == selection_type['all']:
        return slice(None, None, None)
    elif selectionType == selection_type['range']:
        splitString = selectionString.split(":")
        if len(splitString) == 1:
            return slice(int(splitString[0]), None, None)
        elif len(splitString) == 2:
            return slice(int(splitString[0]), int(splitString[1]), None)
        elif len(splitString) == 3: 
            return slice(int(splitString[0]), int(splitString[1]), int(splitString[2]))
        
        
def selection_to_indexlist(selectionString, axisSize=0):
    """Parse the indexlist selection string and return a python list of indicies
    
       :param selectionString: A selection string of the type indexlist
       :param axisSize: Size of the dimensions for which the selection is defined. Only needed in case that a range selection is given.
       
       
       :returns:
            
            * A python list of point indicies for the selection.
            * None in case the list is empty or in case an error occured
    """

    #Check if the given selection is in fact a indexlist
    selectionType = check_selection_string(selectionString) 
    if selectionType == selection_type['indexlist']:
        stringList = selectionString[1:-1].split(",")
        try:
            parsedList = [int(i) for i in stringList]
        except:
            return None
        return parsedList
    elif selectionType == selection_type['index']:
        try:
            return [int(selectionString)]
        except:
            return None
    elif selectionType == selection_type['all']:
        try:
            return range(0, axisSize)
        except:
            return None
    elif selectionType == selection_type['range']:
        try:
            parsedRanges = [int(i) for i in selectionString.split(":")]
            if len(parsedRanges) == 2:
                return range(parsedRanges[0], parsedRanges[1])
            else:
                return None
        except:
            return None
    else:
        return None
    
    


def perform_reduction(data, reduction, axis, httpErr=False):
    """ Helper function used reduce the data of a given numpy array. \n
       
        :param data: The input numpy array that should be reduced
        :param reduction: Data reduction to be applied to the input data. Reduction operations are defined as strings indicating the numpy function to be used for reduction. Valid reduction operations include e.g.: mins, max, mean, median, std, var etc.
        :param axis: The axis along which the reduction should be applied
        :param httpErr: Define which type of error message the function should return. If false then None is returned in case of error. Otherwise a DJANGO HttpResponse is returned.  
        
        :returns: Reduced numpy data array or HttpResonse with a description of the error that occured.
    """
    if httpErr:
        from django.http import HttpResponseNotFound
    
    if data is None:
        return HttpResponseNotFound("Data reduction "+str(reduction)+" failed. None data cannot be reduced.")
    
    axis = int(axis)
    if axis >= len(data.shape):
        return HttpResponseNotFound("Data reduction "+str(reduction)+" failed. The dimensionality of the data is lower than the axis requested to be used for reduction")
    
    #Perform the data reduction operation. This can be a large range of numpy operations eg., min,max,mean,median,std,var
    try:
        op = getattr(np, reduction)
        #if data.shape[axis] > 1:
        data = op(data, axis=axis)
        return data
    except:
        return HttpResponseNotFound("Requested data reduction "+str(reduction)+" failed or not supported. Valid reduction operations are e.g.: min, max, mean, median, std, var."+" "+str(sys.exc_info()))


