import numpy as np
import sys
import itertools


#from omsi.dataformat.omsi_file import omsi_file, omsi_file_analysis
#from omsi.analysis.omsi_analysis_base import omsi_analysis_base

"""This an extended list of types indicated by the check_selection_string function.
   Indicies <0 are assumed to be invalid selections."""
selection_type = {'invalid': -1, 'index': 0, 'indexlist': 2, 'all': 3, 'range': 4}


"""Dicitionary of available data transformation options. Available options are:

   * 'divideMax' : Divide the data by the current maximum value.
   * 'minusMinDivideMax' : Substract the minimum value from the data and \
                   then divide the data by maximum of the data (with the \
                   minimum already substracted.
   * 'logScale' : Apply a log-scale. If the minimum value is 0, then 1 \
                  is added to the data prior to the log scale , ie., np.log(data+1).\
                  If the minimum value is larger then 0, then the log-scale is \
                  applied as is, ie., np.log(data). If the minimum data value is \
                  negative, then the log scale is applied independently to the \
                  positive values and the negative values, ie., \
                  outdata[posvalues] = np.log(data[posvalues]) \
                  outdata[negvalues] = np.log(data[negvalues]*-1.)*-1.
    * 'sqrtScale' : Apply the square-root transformation to the data. The basic \
                  logic in terms of how negative values are treated is the same \
                  as for the logScale transformation.


"""
transformation_type = {'divideMax': 'divideMax',
                       'minusMinDivideMax': 'minusMinDivideMax',
                       'logScale': 'logScale',
                       'sqrtScale': 'sqrtScale'}


def check_selection_string(selection_string):
    """Check whether the given selection string is valid, and indicate which type of selection
       the string defined. Checking the selection string is meant as a safeguard to prevent 
       attakers from being able to insert malitious code.
       
       :param selection_string: String given by the user with the desired selection
       
       :returns: String indicating the type of selection as defined in selection_type:
       
            * 'indexlist' : Selection of the form [1,2,3]
            * 'all' : Selection of the form ':'
            * 'range: Selection of the form 'a:b'
            * 'index: A single index selection, e.g., '1'
            * 'invalid': An unsupported selection 
                 
       """
    import re
    #Check if we actually have a string
    if selection_string is None:
        return selection_type['invalid']
    #Check if the selection defined a list of indicies
    elif re.match('^\[([0-9]+,)*[0-9]+\]$', selection_string):
        return selection_type['indexlist'] 
    #Check if the selection defines and ":" all selection
    elif re.match('^:$', selection_string):
        return selection_type['all']
    #Check if the selection defines a range  a:b type selection
    elif re.match('^[0-9]+:-?[0-9]+$', selection_string):
        return selection_type['range']
    #Check if the selection defines a single index value type selection 
    elif re.match('^-?[0-9]+$', selection_string):
        return selection_type['index']
    #If none of the above selection types are given, then the selection is declared invalid 
    else:
        return selection_type['invalid']


def selection_string_to_object(selection_string):
    """Convert the given selection string to a python selection object, i.e., either a slice, list or integer index.
    
       :param selection_string: A selection string of the type indexlist
       
       :returns: 
       
            * An integer index if an index selection is specified
            * A python list of indices if a list specified in the string
            * A python slice opject if a slice opteration is specified by the string
            
    """
    selectiontype = check_selection_string(selection_string)
    if selectiontype == selection_type['indexlist']:
        stringlist = selection_string[1:-1].split(",")
        try:
            parsedlist = [int(i) for i in stringlist]
        except:
            return None
        return parsedlist
    elif selectiontype == selection_type['index']:
        try: 
            return int(selection_string)
        except:
            return None
    elif selectiontype == selection_type['all']:
        return slice(None, None, None)
    elif selectiontype == selection_type['range']:
        splitstring = selection_string.split(":")
        if len(splitstring) == 1:
            return slice(int(splitstring[0]), None, None)
        elif len(splitstring) == 2:
            return slice(int(splitstring[0]), int(splitstring[1]), None)
        elif len(splitstring) == 3:
            return slice(int(splitstring[0]), int(splitstring[1]), int(splitstring[2]))
        
        
def selection_to_indexlist(selection_string, axis_size=0):
    """Parse the indexlist selection string and return a python list of indices
    
       :param selection_string: A selection string of the type indexlist
       :param axis_size: Size of the dimensions for which the selection is defined.
              Only needed in case that a range selection is given.

       :returns:
            
            * A python list of point indices for the selection.
            * None in case the list is empty or in case an error occurred.

    """

    #Check if the given selection is in fact a indexlist
    selectiontype = check_selection_string(selection_string)
    if selectiontype == selection_type['indexlist']:
        stringlist = selection_string[1:-1].split(",")
        try:
            parsedlist = [int(i) for i in stringlist]
        except:
            return None
        return parsedlist
    elif selectiontype == selection_type['index']:
        try:
            return [int(selection_string)]
        except:
            return None
    elif selectiontype == selection_type['all']:
        try:
            return range(0, axis_size)
        except:
            return None
    elif selectiontype == selection_type['range']:
        try:
            parsedranges = [int(i) for i in selection_string.split(":")]
            if len(parsedranges) == 2:
                return range(parsedranges[0], parsedranges[1])
            else:
                return None
        except:
            return None
    else:
        return None
    

def perform_reduction(data, reduction, axis, http_error=False):
    """ Helper function used reduce the data of a given numpy array.
       
        :param data: The input numpy array that should be reduced
        :param reduction: Data reduction to be applied to the input data.
                          Reduction operations are defined as strings indicating
                          the numpy function to be used for reduction. Valid
                          reduction operations include e.g.: mins, max, mean,
                          median, std, var etc.
        :param axis: The axis along which the reduction should be applied
        :param http_error: Define which type of error message the function should return.
               If false then None is returned in case of error. Otherwise a DJANGO HttpResponse is returned.
        
        :returns: Reduced numpy data array or in case of error None or HttpResonse with a
                  description of the error that occurred (see http_error option).
    """
    if http_error:
        from django.http import HttpResponseNotFound
    
    if data is None:
        if http_error:
            return HttpResponseNotFound("Data reduction "+str(reduction)+" failed. None data cannot be reduced.")
        else:
            return None
    
    axis = int(axis)
    if axis >= len(data.shape):
        if http_error:
            return HttpResponseNotFound("Data reduction "+str(reduction)+" failed." +
                                        "The dimensionality of the data is lower than the axis" +
                                        "requested to be used for reduction")
        else:
            return None
    
    #Perform the data reduction operation. This can be a large range of
    #numpy operations eg., min,max,mean,median,std,var
    try:
        op = getattr(np, reduction)
        #if data.shape[axis] > 1:
        data = op(data, axis=axis)
        return data
    except:
        if http_error:
            return HttpResponseNotFound("Requested data reduction "+str(reduction) +
                                        " failed or not supported. Valid reduction" +
                                        "operations are e.g.: min, max, mean, median, std, var." +
                                        " "+str(sys.exc_info()))
        else:
            return None


def transform_data_multi(data, transformations, http_error=False):
    """ Helper function used to apply a series of potentially multiple
        transformations to a given numpy dataset. This function uses
        the transform_data_single(...) function to apply each indicated
        transformation to the data.

        :param data: The input numpy array that should be transformed.
        :param http_error: Define which type of error message the function should return.
               If false then None is returned in case of error. Otherwise a DJANGO HttpResponse is returned.

        :returns: Reduced numpy data array or HttpResonse with a description of the error that occurred.

    """
    try:
        import json
    except ImportError:
        from django.utils import simplejson as json
    from django.http import HttpResponse, HttpResponseNotFound, HttpResponseBadRequest

    #1) Load the JSON specification of data transformation if necessary
    if isinstance(transformations, str) or isinstance(transformations, unicode):
        try:
            transformations = json.loads(transformations)
        except:
            HttpResponseBadRequest('Invalid list of data transformations.')

    #2) Make sure that we have a list of transformation to iterate over
    if isinstance(transformations, dict):
        transformations = [transformations]

    #3) Iterate over all transformations and apply them on after another.
    for transform in transformations:
        #3.1) Specify the axes for the data transforms
        axes = -1
        if 'axes' in transform:
            axes = transform['axes']
        if 'transformation' in transform:
            currtransformation = transform['transformation']
        else:
            if http_error:
                return HttpResponseBadRequest('Invalid list of transformation. ' +
                                       'Missing transformation key in at least one item')
            else:
                return None
        #3.2) Execute the current data transform
        data = transform_data_single(data=data,
                                     transformation=currtransformation,
                                     axes=axes,
                                     http_error=http_error)
        #3.3) Check if an error occured during the data transformation
        if data is None or isinstance(data, HttpResponse):
            if http_error:
                if isinstance(data, HttpResponse):
                    return data
                else:
                    return HttpResponseNotFound("Unknown error during data transformation.")
            return None

    #4) Return the transformed data
    return data


def transform_data_single(data, transformation=transformation_type['minusMinDivideMax'], axes=-1, http_error=False):
    """ Helper function used to transform data of a numpy array. The function
        potentially splits the array into independent chunks that are
        normalized separately (depending on how the axes parameter is defined).
        The actual data transformations are implemented by transform_datachunk(...).

        :param data: The input numpy array that should be transformed.
        :param transformation: Data transformation option to be used. Available options are:
               'minusMinDivideMax' ,...
        :param axes: List of data axis that should be split into chunks that are treated
                     independently during the transformation. By default transformation is
                     applied based on the full dataset (axes=-1). E.g, if transformation
                     should be performed on a per image basis, then we need to split the
                     m/z dimension into individual chunks and set axes=[2]. If we want
                     to normalize spectra individually, then we need to split the two
                     image dimensions into chunks by setting axes=[0,1].
        :param http_error: Define which type of error message the function should return.
               If false then None is returned in case of error. Otherwise a DJANGO HttpResponse is returned.

        :returns: Reduced numpy data array or HttpResonse with a description of the error that occurred.

    """
    if http_error:
        from django.http import HttpResponseNotFound, HttpResponseBadRequest

    #1) Perform basic error checking
    #1.1) Check input data
    if data is None:
        if http_error:
            return HttpResponseNotFound("Data transformation "+str(transformation) +
                                        " failed. None data cannot be reduced.")
        else:
            return None

    #1.2) Check the axes parameter
    if not isinstance(axes, list):
        axes = [int(axes)]
    for axisindex in axes:
        if axisindex >= len(data.shape):
            if http_error:
                return HttpResponseNotFound("Data transformation "+str(transformation)+" failed." +
                                            " The dimensionality of the data is lower than the axes " +
                                            "requested to be used for reduction")
            else:
                return None

    #1.3) Check the transformation option
    if transformation not in transformation_type:
        if http_error:
            return HttpResponseBadRequest("Data transformation failed. Unsupported transformation option" +
                                          " given as input. Supported options are: "+str(transformation_type))
        else:
            return None

    #2) Normalize the data
    #2.1) Normalize the complete data if no axes is specified
    if axes[0] == -1:
        return transform_datachunk(data=data, transformation=transformation)
    #2.2) Normalize the different chunks based on which axes are specified
    else:
        axislists = []
        for axisindex in axes:
            axislists.append(range(data.shape[axisindex]))
        chunks = itertools.product(*axislists)
        outdata = np.zeros(shape=data.shape, dtype=np.dtype('float'))
        for chunk in chunks:
            selection = [slice(None)]*len(data.shape)
            for coordindex in range(len(axes)):
                selection[axes[coordindex]] = slice(chunk[coordindex],  chunk[coordindex]+1)
            outdata[selection] = transform_datachunk(data=data[selection], transformation=transformation)
        return outdata


def transform_datachunk(data, transformation=transformation_type['minusMinDivideMax']):
    """ Helper function used to transform a given data chunk.
        In contrast to transform_data, this function applies the transformation
        directly to the data provided, without consideration of axis information.
        This function is used by transform_data(...) to implement the actual
        normalization for independent data chunks that need to be normalized.

        :param data: The input numpy array that should be transformed.
        :param transformation: Data transformation option to be used. For available options
                    see the transformation_type dictionary.

        :returns: This function returns the normalized data array. If an unsupported
                  transformation option is given, then the function simply return the
                  unmodified input array.

    """
    if transformation == transformation_type['divideMax']:
        maxvalue = float(np.max(np.abs(data)))
        if maxvalue > 0:
            return data/float(np.max(data))
        else:
            return data
    elif transformation == transformation_type['minusMinDivideMax']:
        minvalue = np.min(data)
        maxvalue = float(np.max(data-minvalue))
        if maxvalue > 0:
            return (data-minvalue) / maxvalue
        else:
            return data
    elif transformation == transformation_type['logScale']:
        minvalue = np.min(data)
        if minvalue == 0:
            return np.log(data+1)
        elif minvalue > 0:
            return np.log(data)
        else:  # minvalue<0
            outdata = np.zeros(shape=data.shape, dtype=np.dtype('float'))
            posvalues = data > 0
            negvalues = data < 0
            outdata[posvalues] = np.log(data[posvalues])
            outdata[negvalues] = np.log(data[negvalues]*-1.)*-1.
            return outdata
    elif transformation == transformation_type['sqrtScale']:
        minvalue = np.min(data)
        if minvalue >= 0:
            return np.sqrt(data)
        else:  # minvalue<0
            outdata = np.zeros(shape=data.shape, dtype=np.dtype('float'))
            posvalues = data >= 0
            negvalues = data < 0
            outdata[posvalues] = np.sqrt(data[posvalues])
            outdata[negvalues] = np.sqrt(data[negvalues]*-1.)*-1.
            return outdata
    else:
        return data


"""
from omsi.shared.omsi_data_selection import *
import json
t = [ {'transformation':'divideMax'} , {'transformation':'logScale', 'axes': [0,1]} ]
tj = json.dumps(t)
a = np.arange(1000).reshape((10,10,10))
at = transform_data_multi(data=a, transformations=tj, http_error=True)
"""