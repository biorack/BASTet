"""
Module for defining and processing data selections. This includes the definition of
selections using strings as well as transformation and reduction of data.

TODO: We may want to expose some of the following numpy functions currently not yet
      supported through the transform and reduce data operations:

    * array2string
    * array_equal
    * array_equiv
    * array_repr
    * array_split
    * array_str
    * asanyarray
    * asarray
    * asarray_chkfinite
    * ascontiguousarray
    * asfarray
    * asfortranarray
    * asmatrix
    * asscalar
    * atleast_1d
    * atleast_2d
    * atleast_3d
    * binary_repr
    * convolve
    * conjugate
    * cross
    * dot
    * extract
    * fft.*
    * histogram, histogram2D, histogramdd
    * kron
    * linalg.*
    * swapaxes(a, axis1, axis2)
    * transpose



#Simple data transformation and reduction example
from omsi.shared.omsi_data_selection import *
import numpy as np
import json
t = [ {'transformation':'threshold', 'threshold':60} , {'reduction':'max', 'axis':2} ]
tj = json.dumps(t)
tj
a = np.arange(125).reshape((5,5,5))
apro = transform_and_reduce_data(data=a, operations=tj, http_error=True)
apro


#Another simple example
from omsi.shared.omsi_data_selection import *
import numpy as np
import json
a = np.arange(10)+5
print a
# 1) substract minimum
# 2) divide by the maximum value with the maximum value converted to float
# NOTE: The conversion to float is to avoid division of integers, i.e.,
#       5/10 = 0, whereas 5/float(10) = 0.5
# NOTE: The specification of 'x1':'data' can be omitted as this is the default.
#       'x1':'data' simply explicitly specifies that the input data should be
#       assigned to the first operand of the arithmetic operation.
t = [{'transformation':'dualDataTransform' , 'operation':'subtract', 'x1':'data', 'x2':[{'reduction':'min'}]} ,
     {'transformation':'dualDataTransform' , 'operation':'divide'  , 'x1':'data', 'x2':[{'reduction':'max'} ,
                                                                    {'transformation':'astype', 'dtype':'float'} ]}]
b = transform_and_reduce_data(data=a, operations=t)
print b
t = [{'transformation':'threshold' , 'threshold':[{'reduction':'median'}]}]
print t
c = transform_and_reduce_data(data=a, operations=t)
print c


#Construct a JSON description of a transformation/reduction
from omsi.shared.omsi_data_selection import *
#Construct the different pieces of the transformation and reduction pipeline
#1) Compute the maximum data value and convert it to float
#1.1) Compute the maximum value
max_value = construct_reduce_dict( reduction_type='max' , axis=None)
#1.2) Convert data to float
value_as_float = construct_transform_dict( trans_type='astype' , dtype='float' )
#1.3) Merge the two steps to compute the maximum data value as float
max_value_as_float = construct_transform_reduce_list( max_value, value_as_float )
#2) Normalize the data by dividing by the maximum value
divide_by_max_value = construct_transform_dict( trans_type='dualDataTransform',
                                                operation='divide' , axes=None , x2=max_value_as_float)
#3) Project along the last axis (i.e., the mz axis) to compute a maximum project image
max_projection = construct_reduce_dict( reduction_type='max' , axis=-1)
#4) Merge the different steps and construct the json string
json_string = transform_reduce_description_to_json( divide_by_max_value , max_projection )
#Just copy the result of the following print statement as your JSON description
print json_string


"""

import numpy as np
import sys
import itertools
import warnings

################################################################
#  Additional custom data processing functions                 #
################################################################
def area_near_peak(x1, halfpeakwidth=2, **kwargs):
    """
    Instead of simply taking the sum or the max to get the "size" of a peak, this code
    finds the peak within a range and sums the data points 'halfpeakwidth' to the
    left and to the right of the peak if there is no peak, the boundary condition ends up being
    that it just takes the first halfpeakwidth+1 datapoints in the range
    but if there is no peak, we assume this value is very low anyways.

    :param x1: The data operand specifying the data the reduction should be performed on
    :param halfpeakwidth: The halfpeakwidth around which the peak should be integrated
    :param kwargs: Additional data parameters. These will be ignored. This is used
                   to support the standard axis parameter that is ignored here.
    """
    if x1.ndim not in [2,3]:
        raise ValueError("Area near peak only defined for 3D cubes of spectra")
    if x1.ndim == 2:
        return x1
    result = np.zeros((x1.shape[0],x1.shape[1]), dtype=x1.dtype)
    maxMasses = np.argmax(x1,2)
    for x in range(x1.shape[0]):
        for y in range(x1.shape[1]):
            peak=(range(-halfpeakwidth, halfpeakwidth+1) + maxMasses[x,y]).astype(int)
            for p in peak:
                if p<0:
                    continue
                if p>= x1.shape[2]:
                    continue
                result[x,y] += x1[x,y,p]
    return result


#################################################################
#  Basic specifications                                         #
#################################################################
reduction_allowed_functions = {'area_near_peak': area_near_peak,
                               'all': np.all,
                               'alltrue': np.alltrue,
                               'amax': np.amax,
                               'amin': np.amin,
                               'angle': np.angle,
                               'any': np.any,
                               'append': np.append,
                               'argmax': np.argmax,
                               'argmin': np.argmin,
                               'argwhere': np.argwhere,
                               'average': np.average,
                               'bincount': np.bincount,
                               'corrcoef': np.corrcoef,
                               'cumprod': np.cumprod,
                               'cumproduct': np.cumproduct,
                               'cumsum': np.cumsum,
                               'count_nonzero': np.count_nonzero,
                               'diag': np.diag,
                               'diag_indices': np.diag_indices,
                               'diagflat': np.diagflat,
                               'diagonal': np.diagonal,
                               'diff': np.diff,
                               'max': np.max,
                               'min': np.min,
                               'median': np.median,
                               'mean': np.mean,
                               'percentile': np.percentile,
                               'product': np.product,
                               'prod': np.prod,
                               'ptp': np.ptp,
                               'squeeze': np.squeeze,
                               'std': np.std,
                               'var': np.var,
                               'transpose': np.transpose,
                               'sum': np.sum,
                               'select_values': None, # Handled separately in the perform_reduction
                               }
"""List of allowed numpy data reduction operations.
   Reduction operations are any single data operations that
   may change the shape of the data.
   NOTE: Some operations may have additional optional or required keyword arguments.
   HELP: For full documentation of the different functions see the numpy documentation.

   Additional input parameters are often:

    * 'x1' : The data operand specifying the data the reduction should be performed on. \
             The input data will be used by default if x1 is not specified. \
             You may also specify 'data' to explicitly indicate that the \
             input data should be assigned to x1. You may specify data0 to \
             indicate that the output of another data operation should be used. \
             Note, data0 here refers to the input to the full data operation \
             pipeline. Data from other parts of the pipeline, are then indexed \
             using 1-based indices. E.g,. to access the output of the first \
             data operation set x1='data0'
    * 'axis' : Integer indicating the axis along which the data should be reduced. \
             The default behavior, if axis is not specified, depends on the behavior \
             of the corresponding numpy function. However, in most cases (if not all cases) \
             the data operation will be applied to the full input data if no axis \
             is specified.
    * 'min_dim' : Integer specifying the minimum number of data dimensions the input data \
             must have in order for the reduction operation to be applied.

   Here the list of allowed data reduction operations.

    * 'all'       : out = numpy.all(data, axis)
    * 'amax'      : out = numpy.amax(data, axis)
    * 'amin'      : out = numpy.amin(data, axis)
    * 'alltrue'   : out = numpy.alltrue(data, axis)
    * 'angle'     : out = numpy.angle(z, deg)
    * 'any'       : out = numpy.any(data, axis)
    * 'append'    : out = numpy.append(data, values, axis)
    * 'argmax'    : out = numpy.argmax(data, axis)
    * 'argmin'    : out = numpy.argmin(data, axis)
    * 'argwhere'  : out = numpy.argwhere(data)
    * 'average'   : out = numpy.average(data, axis)
    * 'bincount'  : out = numpy.bincount(x, weights=None, minlength=None)
    * 'corrcoef'  : out = numpy.corrcoef(data)
    * 'count_nonzero' : out = numpy.count_nonzero(data)
    * 'cumprod'   : out = numpy.cumprod(data,axis)
    * 'cumproduct': out = numpy.cumproduct(data,axis)
    * 'cumsum'    : out = numpy.cumsum(data,axis)
    * 'diag'      : out = numpy.diag(data,k=0)
    * 'diag_indices: out = numpu.diag_indices(data, ndim=2)
    * 'diagflat'  : out = numpy.diagflat(data, k=0)
    * 'diagonal'  : out = numpy.diagonal(data, offset=0, axis1=0, axis2=1)
    * 'diff'      : out = numpy.diff(a, n=1, axis=-1)
    * 'max'       : out = numpy.max(data, axis)
    * 'min'       : out = numpy.min(data, axis)
    * 'median'    : out = numpy.median(data, axis)
    * 'mean'      : out = numpy.mean(data, axis)
    * 'percentile': out = numpy.percentile(data, q, axis)
    * 'product'   : out = numpy.product(data, axis)
    * 'prod'      : out = numpy.prod(data,axis)
    * 'ptp'       : out = numpy.ptp(data,axis)
    * 'squeeze'   : out = numpy.squeeze(data)
    * 'std'       : out = numpy.std(data, axis)
    * 'swapaxes: out = numpy.swapaxes(x1, axis1, axis2)
    * 'var'       : out = numpy.var(data, axis)
    * 'transpose' : out = numpy.transpose(data)
    * 'sum'       : out = numpy.sum(data, axis)

    None-numpy data reduction operations:

    * 'select_values' : out = data[ selection ]

"""

transformation_allowed_single_data = {'abs': np.abs,
                                      'arccos': np.arccos,
                                      'arccosh': np.arccosh,
                                      'arcsin': np.arcsin,
                                      'arcsinh': np.arcsinh,
                                      'arctan': np.arctanh,
                                      'arctanh': np.arctanh,
                                      'argsort': np.argsort,
                                      'around': np.around,
                                      'ceil': np.ceil,
                                      'clip': np.clip,
                                      'cos': np.cos,
                                      'cosh': np.cosh,
                                      'deg2rad': np.deg2rad,
                                      'degrees': np.degrees,
                                      'exp': np.exp,
                                      'exp2': np.exp2,
                                      'fabs': np.fabs,
                                      'floor': np.floor,
                                      'hypot': np.hypot,
                                      'invert': np.invert,
                                      'log': np.log,
                                      'log2': np.log2,
                                      'log10': np.log10,
                                      'logical_not': np.logical_not,
                                      'negative': np.negative,
                                      'sign': np.sign,
                                      'round': np.round,
                                      'sin': np.sin,
                                      'sinc': np.sinc,
                                      'sinh': np.sinh,
                                      'sqrt': np.sqrt,
                                      'sort': np.sort,
                                      'tan': np.tan,
                                      'tanh': np.tanh}
"""List of allowed single data transformations. Single data transformation,
   are operations that operate on a single data input and which do not change
   the shape of the data. Below a list of available numpy options.
   NOTE: Some operations may have additional optional or required keyword arguments.
   HELP: For full documentation of the different functions see the numpy documentation.

   * 'abs'    : out = numpy.abs(x1)
   * 'arccos' : out = numpy.arccos(x1)
   * 'arccosh': out = numpy.arccosh(x1)
   * 'arcsin' : out = numpy.arcsin(x1)
   * 'arcsinh': out = numpy.arcsinh(x1)
   * 'arctan' : out = numpy.arctan(x1)
   * 'arctanh': out = numpy.arctanh1(x1)
   * 'argsort'   : out = numpy.argsort(data, axis, kind='quicksort', order=None)
   * 'around' : out = numpy.around(x1, decimals)
   * 'ceil'   : out = numpy.ceil(x1)
   * 'cos'    : out = numpy.cos(x1)
   * 'cosh'   : out = numpy.cosh(x1)
   * 'clip'   : out = numpy.clip(x1, a_min, a_max)
   * 'deg2rad': out = numpy.deg2rad(x1)
   * 'degrees : out = numpy.degrees(x1)
   * 'exp'    : out = numpy.exp(x1)
   * 'exp2'   : out = numpy.exp2(x1)
   * 'fabs'   : out = numpy.fabs(x1)
   * 'floor'  : out = numpy.floor(x1)
   * 'hypot'  : out = numpy.hypot(x1)
   * 'invert' : out = numpy.invert(x1)
   * 'log'  : out[x1>0] = log(x1[x1>0]) \
              out[x1<0] = log(x1[x1<0]*-1)*-1 \
              out[x1==0] = 0
   * 'log2 ': out[x1>0] = log2(x1[x1>0]) \
              out[x1<0] = log2(x1[x1<0]*-1)*-1 \
              out[x1==0] = 0
   * 'log10': out[x1>0] = log10(x1[x1>0]) \
              out[x1<0] = log10(x1[x1<0]*-1)*-1 \
              out[x1==0] = 0
   * 'logical_not' : out = numpy.logical_not(x1)
   * 'negative' : out = np.negative(x1)
   * 'round' : out = numpy.round(x1, decimals)
   * 'sqrt' : out[x1>0] = sqrt(x1[x1>0]) \
              out[x1<0] = sqrt(x1[x1<0]*-1)*-1 \
              out[x1==0] = 0
   * 'sign' : out = numpy.sign(x1)
   * 'sin'  : out = numpy.sin(x1)
   * 'sinc' : out = numpy.sinc(x1)
   * 'sinh' : out = numpy.sinh(x1)
   * 'sort' : out = numpy.sort(x1, axis=-1, kind='quicksort', order=None)
   * 'swapaxes: out = numpy.swapaxes(x1, axis1, axis2)
   * 'tan'  : out = numpy.tan(x1)
   * 'tanh' : out = numpy.tanh(x1)

"""

transformation_allowed_dual_data = {'add': np.add,
                                    'arctan2': np.arctan2,
                                    'bitwise_and': np.bitwise_and,
                                    'bitwise_not': np.bitwise_not,
                                    'bitwise_or': np.bitwise_or,
                                    'bitwise_xor': np.bitwise_xor,
                                    'corrcoef': np.corrcoef,
                                    'cov': np.cov,
                                    'divide': np.divide,
                                    'equal': np.equal,
                                    'fmax': np.fmax,
                                    'fmin': np.fmin,
                                    'fmod': np.fmod,
                                    'greater': np.greater,
                                    'greater_equal': np.greater_equal,
                                    'left_shift': np.left_shift,
                                    'less': np.less,
                                    'less_equal': np.less_equal,
                                    'logical_and': np.logical_and,
                                    'logical_or': np.logical_or,
                                    'logical_xor': np.logical_xor,
                                    'mod': np.mod,
                                    'multiply': np.multiply,
                                    'not_equal': np.not_equal,
                                    'power': np.power,
                                    'right_shift': np.right_shift,
                                    'subtract': np.subtract}
"""List of allowed dual data transformations. Dual data transformation,
   are operation that operate on a two data input datasets but which do not change
   the shape of the data. Below a list of available numpy function options.
   NOTE: Some operations may have additional optional or required keyword arguments.
   HELP: For full documentation of the different functions see the numpy documentation.

   * 'add'           : out = x1 + x2 = numpy.add(x1,x2)
   * 'arctan2'       : out = numpy.arctan2(x1,x2)
   * 'bitwise_and'   : out = x1 && x2 = numpy.bitwise_and(x1,x2)
   * 'bitwise_not'   : out = numpy.bitwise_not(x1,x2)
   * 'bitwise_or',   : out = x1 || x2 = numpy.bitwise_or(x1,x2)
   * 'bitwise_xor'   : out = numpy.bitwise_xor(x1,x2)
   * 'corrcoef'      : out = numpy.corrcoef(x1,x2)
   * 'cov'           : out = numpy.cov(x1, x2, rowvar=1, bias=0, ddof=None)
   * 'divide'        : out = x1 / x2 = numpy.divide(x1,x2)
   * 'equal'         : out = x1 == x2 = numpy.equal(x1,x2)
   * 'fmax'          : out = numpy.fmax(x1,x2)
   * 'fmin'          : out = numpy.fmin(x1,x2)
   * 'fmod'          : out = numpy.fmod(x1,x2)
   * 'greater'       : out = x1 > x2 = numpy.greater(x1,x2)
   * 'greater_equal' : out = x1 >= x2 = numpy.greater_equal(x1,x2)
   * 'left_shift'    : out = numpy.left_shift(x1,x2)
   * 'less'          : out = x1 < x2 = numpy.less(x1,x2)
   * 'less_equal'    : out = x1 <= x2 = numpy.less_equal(x1,x2)
   * 'logical_and'   : out = numpy.logical_and(x1,x2)
   * 'logical_not'   : See transformation_allowed_single_data instead.
   * 'logical_or'    : out = numpy.logical_or(x1,x2)
   * 'logical_xor'   : out = numpy.logical_xor(x1,x2)
   * 'mod'           : out = numpy.mod(x1,x2)
   * 'multiply'      : out = x1 * x2  = numpy.multiply(x1,x2)
   * 'not_equal'     : out = x1 != x2 = numpy.not_equal(x1,x2)
   * 'power'         : out = numpy.power(x1,x2)
   * 'subtract'      : out = x1 - x2  = numpy.subtract(x1,x2)
   * 'right_shift    : out = np.right_shift(x1,x2)


"""

transformation_type = {'arithmetic': 'arithmetic',
                       'divideMax': 'divideMax',
                       'minusMinDivideMax': 'minusMinDivideMax',
                       'dualDataTransform': 'dualDataTransform',
                       'singleDataTransform': 'singleDataTransform',
                       'scale': 'scale',
                       'threshold': 'threshold',
                       'astype': 'astype'}

"""Dictionary of available data transformation options. Available options are:

   * 'arithmetic' : Same as 'dualDataTransform'. See 'dualDataTransform' below for details.
   * 'divideMax' : Divide the data by the current maximum value.
   * 'minusMinDivideMax' : Substract the minimum value from the data and \
                   then divide the data by maximum of the data (with the \
                   minimum already substracted.
   * 'dualDataTransform' : Apply arbitrary arithmetic operation to the data. Additional parameter \
                  required for this option are:

                    * `operation` : String defining the arithmetic operations to \
                      be applied. Supported operations are: \
                      'add', 'divide', 'greater', 'greater_equal', \
                      'multiply', 'subtract'
                    * 'x1' : The first data operand of the arithmetic operation. \
                             The input data will be used by default if x1 is not specified. \
                             You may also specify 'data' to explicitly indicate that the \
                             input data should be assigned to x1. You may specify data0 to \
                             indicate that the output of another data operation should be used. \
                             Note, data0 here refers to the input to the full data operation \
                             pipeline. Data from other parts of the pipeline, are then indexed \
                             using 1-based indices. E.g,. to access the output of the first \
                             data operation set x1='data0'
                    * 'x2' : The second data operand of the arithmetic operation. \
                             The input data will be used by default if x2 is not specified. \
                             You may also specify 'data' to explicitly indicate that the \
                             input data should be assigned to x2. You may specify data0 to \
                             indicate that the output of another data operation should be used. \
                             Note, data0 here refers to the input to the full data operation \
                             pipeline. Data from other parts of the pipeline, are then indexed \
                             using 1-based indices. E.g,. to access the output of the first \
                             data operation set x2='data0'
                    * ... any additional parameters needed for the numpy function.

   * 'singleDataTransform' : Apply scaling transformation to the data.  Additional parameters \
                  required for this options are. NOTE: operation=='log or operation=='sqrt': \
                  If the minimum value is 0 then the transformation is applied to\
                  positive values only and 0 values remain as is. \
                  If the minimum value is larger then 0, then the log-scale is \
                  applied as is, i.e., np.log(data). If the minimum data value is \
                  negative, then the log scale is applied independently to the \
                  positive values and the negative values, ie., \
                  outdata[posvalues] = np.log(data[posvalues]) \
                  outdata[negvalues] = np.log(data[negvalues]*-1.)*-1.

                    * 'operation' : String defining the scaling operations to be applied. \
                      See the transformation_allowed_single_data list for a complete list
                      of allowed scaling operations. Some of the more commonly used scaling\
                      operations include: 'abs', 'log', 'sqrt', 'around' etc.
                    * 'x1' : The first data operand for the scaling.
                             The input data will be used by default if x1 is not specified. \
                             You may also specify 'data' to explicitly indicate that the \
                             input data should be assigned to x1.

                  Additional optional keyword arguments depending on the used operation:

                    * 'decimals' : Number of decimal places to round to when using numpy.around or numpy.round \
                                  (default: 0).  If decimals is negative, it specifies the
                                  number of positions to the left of the decimal point.
                    * 'a_min', 'a_max' : Lower and upper bound when using numpy.clip.
                    * 'axis', 'kind', 'order' : Additional optional arguments for numpy.argsort and numpy.sort.
                    * ...

    * 'scale' : Same as 'singleDataTransform'. See 'singleDataTransform' for details.
    * 'threshold' : Threshold the data. Set all values that are smaller than threshold \
                  to 0. Additional parameters required for this option are:

                    * 'threshold'. If threshold is missing, then the threshold will be \
                     set ot the 5%'ile so that the bottom 5% of the data will be set to 0.

    * 'astype' : Change the type of the data. Additional required parameters are: \

                    * 'dtype' : The numpy data type to be used. Default dtype='float'.


"""


selection_type = {'invalid': -1,
                  'index': 0,
                  'indexlist': 2,
                  'all': 3,
                  'range': 4,
                  'multiaxis': 5}
"""This an extended list of types indicated by the check_selection_string function.
   Indices <0 are assumed to be invalid selections.

"""


#################################################################
#  Process data selection strings                               #
#################################################################
def check_selection_string(selection_string):
    """
    Check whether the given selection string is valid, and indicate which type of selection
    the string defined. Checking the selection string is meant as a safeguard to prevent
    attackers from being able to insert malicious code.

    :param selection_string: String given by the user with the desired selection

    :returns: Interger indicating the type of selection as defined in selection_type, e.g,:

        * 2 = 'indexlist' : Selection of the form [1,2,3]
        * 3 = 'all' : Selection of the form ':'
        * 4 = 'range: Selection of the form 'a:b'
        * 0 = 'index: A single index selection, e.g., '1'
        * 5 = 'multiaxis': A selection for multiple axies, e.g., '1:2|3:4'
        * -1 'invalid': An unsupported selection

    """
    import re
    # Check if we actually have a string
    if selection_string is None:
        return selection_type['invalid']
    # Check if the selection defined a list of indicies
    elif re.match('^\[([0-9]+,)*[0-9]+\]$', selection_string):
        return selection_type['indexlist']
    # Check if the selection defines and ":" all selection
    elif re.match('^:$', selection_string):
        return selection_type['all']
    # Check if the selection defines a range/slice  a:b:c type selection
    elif re.match('^[0-9]+:-?[0-9]+$', selection_string):
        return selection_type['range']
    # Check if the selection defines a single index value type selection
    elif re.match('^-?[0-9]+$', selection_string):
        return selection_type['index']
    # Check if we have a multi-axis selection
    elif "|" in selection_string:
        for substring in selection_string.split('|'):
            if check_selection_string(substring) < 0:
                return selection_type['invalid']
        return selection_type['multiaxis']
    # Check if we have a more complex slice selection
    elif ":" in selection_string:
        splitstring = selection_string.split(":")
        try:
            slice_values = [int(i) if len(i) > 0 else None for i in splitstring]
        except ValueError:
            return selection_type['invalid']
        contains_digit = False
        for current_slice_value in slice_values:
            if current_slice_value is not None:
                contains_digit = True
                break
        if (0 < len(slice_values) < 4) and contains_digit:
            return selection_type['range']
        else:
            return selection_type['invalid']
    # If none of the above selection types are given, then the selection is
    # declared invalid
    else:
        return selection_type['invalid']


def selection_string_to_object(selection_string,
                               list_to_index=False):
    """
    Convert the given selection string to a python selection object, i.e., either a slice, list or integer index.

    :param selection_string: A selection string of the type indexlist
    :param list_to_index: Should we turn the list into an index if the list contains only a single value.
        Default value is False, i.e., the list is not modified.

    :returns:

        * An integer index if an index selection is specified
        * A python list of indices if a list specified in the string
        * A python slice object if a slice operation is specified by the string

    """
    selectiontype = check_selection_string(selection_string)
    if selectiontype == selection_type['indexlist']:
        stringlist = selection_string[1:-1].split(",")
        try:
            parsedlist = [int(i) for i in stringlist]
        except ValueError:
            return None
        if len(parsedlist) == 1 and list_to_index:
            return parsedlist[0]
        else:
            return parsedlist
    elif selectiontype == selection_type['index']:
        try:
            return int(selection_string)
        except ValueError:
            return None
    elif selectiontype == selection_type['all']:
        return slice(None, None, None)
    elif selectiontype == selection_type['range']:
        splitstring = selection_string.split(":")
        slice_values = [int(i) if len(i) > 0 else None for i in splitstring]

        if len(splitstring) == 1:
            return slice(slice_values[0], None, None)
        elif len(splitstring) == 2:
            return slice(slice_values[0], slice_values[1], None)
        elif len(splitstring) == 3:
            return slice(slice_values[0], slice_values[1], slice_values[2])
        else:
            return slice(None, None, None)
    elif selectiontype == selection_type['multiaxis']:
        splitstring = selection_string.split('|')
        return tuple([selection_string_to_object(axis_selection) for axis_selection in splitstring])
    elif isinstance(selection_string, list) or \
            isinstance(selection_string, int) or \
            isinstance(selection_string, tuple) or \
            isinstance(selection_string, slice):
        return selection_string
    else:
        return None


def selection_to_string(selection):
    """
    Convert the given selection, which may be either an int, a list of ints, a slice object or
    a tuple of the mentioned types which is used to define a selection along multiple axes.
    :param selection: The selection to be converted to a string
    :type selection: int, list, slice, or a tuple of int, list, slice objects
    :return: The selection string
    """
    if isinstance(selection, list) or isinstance(selection, int):
        return str(selection).replace(' ', '')
    elif isinstance(selection, slice):
        start = selection.start
        stop = selection.stop
        step = selection.step
        if start is None and stop is None and step is None:
            return ":"
        else:
            start_string = "" if start is None else str(start)
            stop_string = "" if stop is None else str(stop)
            step_string = "" if step is None else str(step)
            return "%s:%s:%s" % (start_string, stop_string, step_string)
    elif isinstance(selection, tuple):
        selection_string = None
        for axis_selection in selection:
            if selection_string is None:
                selection_string = selection_to_string(axis_selection)
            else:
                selection_string += "|" + selection_to_string(axis_selection)
        if selection_string is None:
            selection_string = ""
        return selection_string
    elif isinstance(selection, basestring):
        if check_selection_string(selection) >= 0:
            return selection
        else:
            raise ValueError('Invalid selection string given to selection_to_string(...)')
    else:
        raise ValueError('Invalid selection given to selection_to_string(...).')


def selection_to_indexlist(selection_string, axis_size=0):
    """
    Parse the indexlist selection string and return a python list of indices

    :param selection_string: A selection string of the type indexlist
    :param axis_size: Size of the dimensions for which the selection is defined.
          Only needed in case that a range selection is given. This should be
          a list of sizes, in case that a multiaxis selection is given.

    :returns:

        * A python list of point indices for the selection.
        * None in case the list is empty or in case an error occurred.
    """
    # Check if the given selection is in fact a indexlist
    selectiontype = check_selection_string(selection_string)
    if selectiontype == selection_type['indexlist']:
        return selection_string_to_object(selection_string)
    elif selectiontype == selection_type['index']:
        return [selection_string_to_object(selection_string)]
    elif selectiontype == selection_type['all']:
        return range(0, axis_size)
    elif selectiontype == selection_type['range']:
        slice_object = selection_string_to_object(selection_string)
        if slice_object.start is not None:
            start = slice_object.start
        else:
            start = 0
        if slice_object.stop is not None:
            stop = slice_object.stop
        else:
            stop = axis_size
        if slice_object.step is not None:
            step = slice_object.step
        else:
            step = 1
        return range(start, stop, step)
    elif selection_type == selection_type['multiaxis']:
        selection_tuple = selection_string_to_object(selection_string)
        selection_lists = []
        if axis_size is None:
            axis_sizes = [0] * len(selection_tuple)
        elif isinstance(axis_size, list) or isinstance(axis_size, tuple):
            axis_sizes = list(axis_size)
        else:
            axis_sizes = [axis_size] * len(selection_tuple)
        for axis_index, axis_selection in enumerate(selection_tuple):
            selection_lists.append(selection_to_indexlist(selection_string=axis_selection,
                                                          axis_size=axis_sizes[axis_index]))
            return zip(*selection_lists)
    else:
        return None


#################################################################
#  Evaluate data reductions and transformations                #
#################################################################
def perform_reduction(data, reduction, secondary_data, min_dim=None, http_error=False, **kwargs):
    """ Helper function used reduce the data of a given numpy array.

        :param data: The input numpy array that should be reduced
        :param reduction: Data reduction to be applied to the input data.
                          Reduction operations are defined as strings indicating
                          the numpy function to be used for reduction. Valid
                          reduction operations include e.g.: mins, max, mean,
                          median, std, var etc.
        :type reduction: String
        :param axis: The axis along which the reduction should be applied
        :param secondary_data: Other data from previous data iterations a user may reference.
        :param http_error: Define which type of error message the function should return.
               If false then None is returned in case of error. Otherwise a DJANGO HttpResponse is returned.
        :param min_dim: Minimum number of dimensions the input data must have in order for the reduction to be applied.
        :param kwargs: Additional optional keyword arguments.

        :returns: Reduced numpy data array or in case of error None or HttpResonse with a
                  description of the error that occurred (see http_error option).
    """
    try:
        from django.http import HttpResponseNotFound
    except ImportError:
        http_error = False

    # 1) Check input
    # 1.1) Check if we have any data that can be reduced
    if data is None:
        return None

    # 1.2) Check if we have a valid reduction operation
    if reduction not in reduction_allowed_functions:
        if http_error:
            return HttpResponseNotFound("Requested data reduction " + str(reduction) +
                                        " not supported. Valid reduction operations are:" +
                                        str(reduction_allowed_functions.keys()))
        else:
            raise ValueError("Requested data reduction " + str(reduction) +
                             " not supported. Valid reduction operations are:" +
                             str(reduction_allowed_functions.keys()))

    # 1.3) Check if we have a valid axis parameter
    axis_specified = 'axis' in kwargs
    if axis_specified:
        axis = kwargs.pop('axis')
        if axis is not None:
            axis = int(axis)
        if axis >= len(data.shape):
            if http_error:
                return HttpResponseNotFound("Data reduction " + str(reduction) + " failed." +
                                            "The dimensionality of the data is lower than the axis" +
                                            "requested to be used for reduction")
            else:
                raise IndexError("Data reduction " + str(reduction) + " failed." +
                                 "The dimensionality of the data is lower than the axis" +
                                 "requested to be used for reduction")

    # 1.4) Determine the data we need to perform the reduction for
    if 'x1' in kwargs:
        x1 = evaluate_transform_parameter(parameter=kwargs.pop('x1'),
                                          data=data,
                                          secondary_data=secondary_data)
    else:
        x1 = data

    # 1.5) Determine whether we should apply the reduction or not
    if min_dim:
        if x1.ndim < min_dim:
            return x1

    # 1.6) Check if we have a data selection reduction
    if reduction == 'select_values':
        if 'selection' in kwargs:
            selection = kwargs.pop('selection')
            if is_transform_or_reduce(selection):
                selection = evaluate_transform_parameter(parameter=selection,
                                                         data=x1,
                                                         secondary_data=secondary_data)
            if isinstance(selection, str) or isinstance(selection, unicode):
                selection = selection_string_to_object(selection)
        else:
            selection = slice(None)
        return data[selection]

    # 1.7 ) Perform the data reduction operation. This can be a large range of
    # numpy operations defined in
    # omsi.shared.omsi_data_selection.reduction_allowed_functions
    try:
        reduction_operation = reduction_allowed_functions[reduction]  # getattr(np, reduction)
        # if data.shape[axis] > 1:
        if axis_specified:
            data = reduction_operation(x1, axis=axis, **kwargs)
        else:
            data = reduction_operation(x1, **kwargs)
        return data
    except:
        if http_error:
            return HttpResponseNotFound("Requested data reduction " + str(reduction) +
                                        " failed or not supported. Valid reduction" +
                                        "operations are:" + str(reduction_allowed_functions.keys()) +
                                        " " + str(sys.exc_info()))
        else:
            raise ValueError("Requested data reduction " + str(reduction) +
                             " failed or not supported. Valid reduction" +
                             "operations are:" + str(reduction_allowed_functions.keys()) +
                             " " + str(sys.exc_info()))


def transform_and_reduce_data(data,
                              operations,
                              secondary_data=None,
                              http_error=False):
    """ Helper function used to apply a series of potentially multiple
        operations to a given numpy dataset. This function uses
        the transform_data_single(...) function to apply each indicated
        transformation to the data. This function uses the perform_reduction function to
        perform data reduction operations.

        :param data: The input numpy array that should be transformed.
        :param operations: JSON string with list of dictionaries or a python
                           list of dictionaries. Each dict specifies a single data
                           transformation or data reduction. The operations are applied
                           in order, i.e., operations[0] is applied first, then operations[1]
                           and so on. The dicts must be structured according to one of
                           the following specifications:

                           * `{'transformation':<op>}` : Single transformation applied to all data at once.
                           * `{'transformation':<op>, 'axes':[..]}` : Apply a single transformation to \
                             data chunks defined by the axes parameter. The data is split into chunks along \
                             the dimensions defined by the axes parameter. E.g., if we have a 3D MSI dataset \
                             and we want to op ion images independently, then we need to set axes=[2]. \
                             Accordingly,  if we want to op spectra individually, then we need to split \
                             the two image dimensions into chunks by setting axes=[0,1].
                           * `{'reduction':<reduction>, 'axis':int}` : Define the reduction operations to be \
                             applied and the axis along which the data should be reduced. If reduction along \
                             all axis should be done then set axis ot None (in python) or null in JSON.

        :param secondary_data: Other data from previous data iterations a user may reference.
        :param http_error: Define which type of error message the function should return.
               If false then None is returned in case of error. Otherwise a DJANGO HttpResponse is returned.

        :returns: Reduced numpy data array or HttpResonse with a description of the error that occurred.

    """
    try:
        import json
    except ImportError:
        from django.utils import simplejson as json
    from django.http import HttpResponse, HttpResponseNotFound, HttpResponseBadRequest

    # 1) Check if there is anything to do. If not, then return the data.
    if operations is None:
        return data
    if isinstance(operations, list):
        if len(operations) == 0:
            return data

    # 2) Load the JSON specification of data transformation if necessary
    if isinstance(operations, str) or isinstance(operations, unicode):
        try:
            operations = json.loads(operations)
        except:
            HttpResponseBadRequest('Invalid list of data operations.')

    # 3) Make sure that we have a list of transformation to iterate over
    if isinstance(operations, dict):
        operations = [operations]

    # 4) Determine which data outputs we need to track, i.e., check which
    #    data outputs are explicitly referenced
    if secondary_data is None:
        secondary_data = {'data': data}
    elif secondary_data == {}:
        secondary_data['data'] = data

    # 4) Iterate over all operations and apply them on after another.
    current_data = data

    for current_operation in operations:
        # print "---------------------------"
        # print op
        if 'variable' in current_operation:
            secondaryvar = current_operation.pop('variable')
        else:
            secondaryvar = None
        # 4.1) Evaluate a data transformation that consists of multiple
        # operations
        if isinstance(current_operation, list):
            transform_and_reduce_data(data=current_data,
                                      operations=current_operation,
                                      secondary_data=secondary_data,
                                      http_error=http_error)
        # 4.2) Evaluate single data operation
        elif isinstance(current_operation, dict):
            # 4.2.1) Check if we need to transform or reduce the data
            # 4.2.2) Perform data transformation
            if 'transformation' in current_operation:
                # 4.2.2.1 Get the transformation
                current_transformation = current_operation.pop('transformation')
                # 4.2.2.2 Get the axes for the transformations
                axes = None
                if 'axes' in current_operation:
                    axes = current_operation.pop('axes')
                # 4.2.3) Execute the current data operation
                current_data = transform_data_single(data=current_data,
                                                     transformation=current_transformation,
                                                     axes=axes,
                                                     secondary_data=secondary_data,
                                                     http_error=http_error,
                                                     transform_kwargs=current_operation)
            # 4.2.3 Perform data reduction
            elif 'reduction' in current_operation:
                # 4.2.3.1 Get the data reduction operation
                currreduction = current_operation.pop('reduction')
                # 4.2.3.2 Get the axis along which the data reduction should be
                # applied
                if 'axis' in current_operation:
                    axis = current_operation.pop('axis')
                else:
                    axis = None  # data.ndim-1
                if 'min_dim' in current_operation:
                    min_dim = current_operation.pop('min_dim')
                else:
                    min_dim = None
                # 4.2.3.3 Execute the data reduction operation
                current_data = perform_reduction(data=current_data,
                                                 reduction=currreduction,
                                                 secondary_data=secondary_data,
                                                 axis=axis,
                                                 min_dim=min_dim,
                                                 http_error=http_error,
                                                 **current_operation)
            # 4.2.4 Invalid operation specification. We have neither a
            # transformation nor a reduction
            else:
                if http_error:
                    return HttpResponseBadRequest('Invalid list of transformation/reduction operations. ' +
                                                  'Missing transformation or reduction key in at least one item')
                else:
                    raise KeyError('Invalid list of transformation/reduction operations. ' +
                                   'Missing transformation or reduction key in at least one item')

        # 4.3) Check if an error occurred during the data transformation or
        # reduction
        if current_data is None or isinstance(current_data, HttpResponse):
            if http_error:
                if isinstance(current_data, HttpResponse):
                    return current_data
                else:
                    return HttpResponseNotFound("Unknown error during data transformation.")
            return current_data

        # 4.4) Track the secondary data if needed
        if secondaryvar:
            secondary_data[secondaryvar] = current_data

        # print "---------------------------"

    # 5) Return the transformed data
    return current_data


def transform_data_single(data,
                          transformation=transformation_type['minusMinDivideMax'],
                          axes=None,
                          secondary_data=None,
                          http_error=False,
                          transform_kwargs=None):
    """ Helper function used to transform data of a numpy array. The function
        potentially splits the array into independent chunks that are
        normalized separately (depending on how the axes parameter is defined).
        The actual data transformations are implemented by transform_datachunk(...).

        :param data: The input numpy array that should be transformed.
        :param transformation: Data transformation option to be used. Available options are:
               'minusMinDivideMax' ,...
        :param axes: List of data axis that should be split into chunks that are treated
                     independently during the transformation. By default transformation is
                     applied based on the full dataset (axes=None). E.g, if transformation
                     should be performed on a per image basis, then we need to split the
                     m/z dimension into individual chunks and set axes=[2]. If we want
                     to transform spectra individually, then we need to split the two
                     image dimensions into chunks by setting axes=[0,1].
        :param secondary_data: Other data from previous data iterations a user may reference.
        :param http_error: Define which type of error message the function should return.
               If false then None is returned in case of error. Otherwise a DJANGO HttpResponse is returned.
        :param transform_kwargs: Dictionary of additional keyword arguments to be passed to the
                       transform_datachunk(...) function.

        :returns: Reduced numpy data array or HttpResonse with a description of the error that occurred.

    """
    if http_error:
        from django.http import HttpResponseNotFound, HttpResponseBadRequest

    # 1) Perform basic error checking
    # 1.1) Check input data
    if data is None:
        if http_error:
            return HttpResponseNotFound("Data transformation " + str(transformation) +
                                        " failed. None data cannot be reduced.")
        else:
            raise ValueError("Data transformation " + str(transformation) +
                             " failed. None data cannot be reduced.")

    # 1.2) Check the axes parameter
    if not isinstance(axes, list) and axes is not None:
        axes = [int(axes)]
    if axes is not None:
        for axisindex in axes:
            if axisindex >= len(data.shape):
                if http_error:
                    return HttpResponseNotFound("Data transformation " + str(transformation) + " failed." +
                                                " The dimensionality of the data is lower than the axes " +
                                                "requested to be used for reduction")
                else:
                    raise IndexError("Data transformation " + str(transformation) + " failed." +
                                     " The dimensionality of the data is lower than the axes " +
                                     "requested to be used for reduction")

    # 1.3) Check the transformation option
    if transformation not in transformation_type:
        if http_error:
            return HttpResponseBadRequest("Data transformation failed. Unsupported transformation option" +
                                          " given as input. Supported options are: " + str(transformation_type))
        else:
            raise ValueError("Data transformation failed. Unsupported transformation option" +
                             " given as input. Supported options are: " + str(transformation_type))

    # 1.4) Check the additional keyword arguments
    if transform_kwargs is None:
        transform_kwargs = {}

    # 2) Transform the data
    # 2.1) Transform the complete data if no axes is specified
    if axes is None:
        return transform_datachunk(data=data,
                                   transformation=transformation,
                                   secondary_data=secondary_data,
                                   **transform_kwargs)

    # 2.2) Transform the different chunks based on which axes are specified
    else:
        axislists = []
        for axisindex in axes:
            axislists.append(range(data.shape[axisindex]))
        chunks = itertools.product(*axislists)
        outdata = np.zeros(shape=data.shape, dtype=np.dtype('float'))
        for chunk in chunks:
            selection = [slice(None)] * len(data.shape)
            for coordindex in range(len(axes)):
                selection[axes[coordindex]] = slice(
                    chunk[coordindex], chunk[coordindex] + 1)
            outdata[selection] = transform_datachunk(data=data[selection],
                                                     transformation=transformation,
                                                     secondary_data=secondary_data,
                                                     **transform_kwargs)
        return outdata


def transform_datachunk(data,
                        transformation=transformation_type['minusMinDivideMax'],
                        secondary_data=None,
                        **kwargs):
    """ Helper function used to transform a given data chunk.
        In contrast to transform_data, this function applies the transformation
        directly to the data provided, without consideration of axis information.
        This function is used by transform_data(...) to implement the actual
        normalization for independent data chunks that need to be normalized.

        Required keyword arguments:

        :param data: The input numpy array that should be transformed.
        :param transformation: Data transformation option to be used. For available options
                    see the transformation_type dictionary.
        :param secondary_data: Other data from previous data iterations a user may reference.

        Additional transformation-dependent keyword arguments:

        :param kwargs: Additional keyword arguments that are specific for different
                       data transformation. Below a list of additional keyword arguments
                       used for different transformation options

                       * transformation: 'threshold'
                            ** 'threshold' : The threshold parameter to be used for
                                             the thresold operation. If threshold is
                                             not specified, then the 5th %tile will
                                             be used as threshold value instead, ie.,
                                             the bottom 5% of the data are set to 0.

        :returns: This function returns the normalized data array. If an unsupported
                  transformation option is given, then the function simply return the
                  unmodified input array.

    """
    ############################################
    #   Divide by max                          #
    ############################################
    if transformation == transformation_type['divideMax']:
        maxvalue = float(np.max(np.abs(data)))
        if maxvalue > 0:
            return data / float(np.max(data))
        else:
            return data
    ############################################
    #   minus min divide by max                #
    ############################################
    elif transformation == transformation_type['minusMinDivideMax']:
        minvalue = np.min(data)
        maxvalue = float(np.max(data - minvalue))
        if maxvalue > 0:
            return (data - minvalue) / maxvalue
        else:
            return data
    ############################################
    #   arithmetic                             #
    ############################################
    elif transformation == transformation_type['dualDataTransform'] or \
            transformation == transformation_type['arithmetic']:
        # Retrieve the requested operation
        if 'operation' in kwargs:
            if kwargs['operation'] not in transformation_allowed_dual_data:
                raise AttributeError('Requested data operations not supported.')
            try:
                operation = transformation_allowed_dual_data[kwargs.pop('operation')] # getattr(np, kwargs.pop('operation'))
            except AttributeError:
                raise AttributeError('Requested data operations not supported.')
        else:
            raise KeyError("Missing operation key for arithmetic data transformation.")
        # Evalute the parameters x1 and x2
        if 'x1' in kwargs:
            x1 = evaluate_transform_parameter(parameter=kwargs.pop('x1'),
                                              data=data,
                                              secondary_data=secondary_data)
        else:
            x1 = data
        if 'x2' in kwargs:
            x2 = evaluate_transform_parameter(parameter=kwargs.pop('x2'),
                                              data=data,
                                              secondary_data=secondary_data)
        else:
            x2 = data
        # Execute the requested data operation
        return operation(x1, x2, **kwargs)
    ############################################
    #   scale                                  #
    ############################################
    elif transformation == transformation_type['singleDataTransform'] or transformation == transformation_type['scale']:
        specialzero = False
        if 'operation' in kwargs:
            if kwargs['operation'] not in transformation_allowed_single_data:
                raise AttributeError(
                    'Requested scaling operation not supported.')
            try:
                if kwargs['operation'] == 'log' or \
                   kwargs['operation'] == 'log10' or \
                   kwargs['operation'] == 'sqrt':
                    specialzero = True
                operation =  transformation_allowed_single_data[kwargs.pop('operation')]  # getattr(np, kwargs.pop('operation'))
            except AttributeError:
                raise AttributeError(
                    'Requested data operations not supported.')
        else:
            raise KeyError(
                "Missing operation key for scaling data transformation.")
        # Evalute the parameters x1 and x2
        if 'x1' in kwargs:
            x1 = evaluate_transform_parameter(parameter=kwargs.pop('x1'),
                                              data=data,
                                              secondary_data=secondary_data)
        else:
            x1 = data
        # Execute the scaling operation operation
        if not specialzero:
            return operation(x1, **kwargs)
        else:
            minvalue = np.min(data)
            if minvalue > 0:
                return operation(x1, **kwargs)
            elif minvalue <= 0:
                outdata = np.zeros(shape=data.shape, dtype=np.dtype('float'))
                posvalues = x1 > 0
                negvalues = x1 < 0
                outdata[posvalues] = operation(x1[posvalues], **kwargs)
                outdata[negvalues] = operation(
                    x1[negvalues] * -1, **kwargs) * -1.
                return outdata
    ############################################
    #   Threshold                              #
    ############################################
    elif transformation == transformation_type['threshold']:
        outdata = np.copy(data)
        if 'threshold' in kwargs:
            threshold = evaluate_transform_parameter(
                parameter=kwargs.pop('threshold'),
                data=data,
                secondary_data=secondary_data)
        else:
            threshold = np.percentile(data, 0.05)
        outdata[data < threshold] = 0
        return outdata
    ############################################
    #   astype                                 #
    ############################################
    elif transformation == transformation_type['astype']:
        if 'dtype' in kwargs:
            dtype = np.dtype(kwargs.pop('dtype'))
        else:
            dtype = np.dtype('float')
        if 'x1' in kwargs:
            x1 = evaluate_transform_parameter(parameter=kwargs.pop('x1'),
                                              data=data,
                                              secondary_data=secondary_data)
        else:
            x1 = data
        return x1.astype(dtype)
    else:
        return data


def evaluate_transform_parameter(parameter,
                                 data=None,
                                 secondary_data=None):
    """Evaluate the given query parameter. This function is used to
       enable the use of data transformation and reductions as part
       of transformation parameters. E.g., a user may want to substract
       the minimum, or divide by the maximum etc.

       :param parameter: The parameter to be evaluated. This may be
              a JSON string or list/dictionary-based description of a
              data transformation. Or any other valid data parameter.
              If the parameter describes as data reduction or transformation
              then the transformation will be evaluated and the result
              is returned, otherwise the parameter itself is returned.
       :param data: The input numpy array that should be transformed.
       :param secondary_data: Other data from previous data iterations a user may reference.

       :returns: The evaluated parameter result.
    """
    if is_transform_or_reduce(parameter):
        return transform_and_reduce_data(data=data,
                                         operations=parameter,
                                         secondary_data=secondary_data,
                                         http_error=False)
    else:
        if parameter in secondary_data or unicode(parameter) in secondary_data:
            return secondary_data[parameter]
        else:
            return parameter


#################################################################
#   Helper functions to assist with the construction and        #
#   validation of data operation descriptions for data          #
#   transformation and reduction.                               #
#################################################################
# def get_secondary_data_used( operations , optimize_specification=False ):
#     """Check for the all operations in the list if they directly reference
#        any data from any of the other steps. This check is only for the
#        main level operations, i.e., those operations in the list itself and
#        not for sub-operations that may be part of any of the parameters of
#        data operations. This function is used to mark which data outputs
#        should be tracked as part of the computation.
#
#        :param operations: List of operations to be performed.
#        :param optimize_specification: If set to True then the function tries to
#               optimize the given specification of operations by removing any
#               unneeded secondary data dependencies. E.g., if any operation
#               explicitly references the data from the previous step, then
#               the use of secondary data can be avoided by replacing the
#               reference with a simple reference to the data itself. This
#               helps to reduce the amount of secondary data that needs to be
#               tracked.
#
#        :raises: IndexError in case a secondary data-index is invalid
#
#     """
#     from copy import copy
#     if isinstance(operations,dict):
#         oplist = [copy(operations)]
#     else:
#         oplist = copy(operations)
#
#     secondaryused = np.zeros( len(operations)+1 , dtype=np.dtype('bool') )
#     index = 1
#     for op in oplist:
#         if isinstance(op,dict):
#             for key, value in op.items():
#                 if isinstance(value,str) or isinstance(value,unicode):
#                     if value.startswith('data'):
#                         try:
#                             reference = int(value.lstrip('data'))
#                             if reference == (index-1) and optimize_specification:
#                                 op[key] = 'data'
#                             else:
#                                 secondaryused[reference] = True
#                         except ValueError:
#                             pass
#                         except IndexError:
#                             raise
#         index += 1
#     return secondaryused, oplist


def is_transform_or_reduce(parameter):
    """Check if the given parameter defines a description of a
       data transformation or data reduction

       :param parameter: The parameter to be checked.
       :type parameter: JSON string, dict or list of dicts with transformation parameter.

       :returns: Boolean
    """
    try:
        import json
    except ImportError:
        from django.utils import simplejson as json

    # 1) If the parameter is a built in type, then return the parameter
    if isinstance(parameter, int) or \
       isinstance(parameter, float) or \
       isinstance(parameter, long) or \
       isinstance(parameter, complex) or \
       isinstance(parameter, bool):
        return False

    # 2) Try to convert the parameter to the list of dict based description
    evalparam = parameter
    if isinstance(parameter, str) or isinstance(parameter, unicode):
        if parameter == 'data':
            return False
        try:
            evalparam = json.loads(parameter)
        except ValueError:
            pass
        except TypeError:
            pass
    if isinstance(evalparam, dict):
        evalparam = [evalparam]

    # 3) Check if we have a valid description
    if isinstance(evalparam, list):
        # Check if each element has the minimum set of parameters in the dicts
        for paramelement in evalparam:
            if isinstance(paramelement, dict):
                if 'transformation' in paramelement or \
                   'reduction' in paramelement:
                    pass
                else:
                    return False
            else:
                return False
            return True
    else:
        return False


def construct_transform_reduce_list(*args):
    """Merge a series of transformations and reductions into a single
       list describing a pipeline of transformation and reduction
       operations to be performed.

       :args: Ordered series of dictionaries describing transformation
              and reduction operations.

       :returns: List of all transformation and reduction operations
    """
    # Check if all elements are valid
    count = 0
    for element in args:
        if not is_transform_or_reduce(element):
            warnings.warn("Transformation/Reduction " + str(count) +
                          " may be invalid: " + str(element))
        count += 1
    # Convert the arguments to a list and return
    return list(args)


def construct_transform_dict(trans_type, axes=None, **kwargs):
    """Helper function used to construct a dictionary describing a data transformation.

       :param trans_type: The transformation type to be used. See transformation_type dict.
       :param axes: The axes along which the data should be split. Default is None.
       :param kwargs: Additional keyword parameters for the transformation functions.

       :returns: Dictionary with the description of the transformation.

       :raises: KeyError is raised in case that a parameter is missing.
                ValueError is raised in case that a given parameter value
                is invalid.

    """
    # 1) Validate input parameters
    # 1.1) Do we have a valid transformation type
    if trans_type not in transformation_type:
        raise ValueError('Invalid transformation type given. Valid trans_type are:' +
                         str(transformation_type))
    # 1.1) Do we have a valid operation specified for arithmetic opertions
    if trans_type == transformation_type['dualDataTransform'] or trans_type == transformation_type['arithmetic']:
        if 'operation' in kwargs:
            if not kwargs['operation'] in transformation_allowed_dual_data:
                raise ValueError(unicode("Requested arithmetic operation not supported. Allowed operations are") +
                                 unicode(transformation_allowed_dual_data.keys()))
        else:
            raise KeyError(
                "Missing parameter operations for arithmetic data transformation.")
    # 1.2) Check if we have a scaling operation and if so if it is valid
    if trans_type == transformation_type['singleDataTransform'] or trans_type == transformation_type['scale']:
        if 'operation' in kwargs:
            if not kwargs['operation'] in transformation_allowed_single_data:
                raise ValueError(unicode("Requested scale operation not supported. Allowed operations are") +
                                 unicode(transformation_allowed_single_data.keys()))
        else:
            raise KeyError("Missing parameter operations for scale data transformation.")
    # 1.3) Check if we have threshold parameter for the threshold operation
    if trans_type == transformation_type['threshold']:
        if 'threshold' not in kwargs:
            warnings.warn('No threshold parameter specified. The 5th percentile will be used as default')
    # 1.4) Check whether we have a dtype for the astype function
    if trans_type == transformation_type['astype']:
        if 'dtype' not in kwargs:
            raise KeyError("Missing parameter dtype for astype transformation")

    transdict = {'transformation': unicode(trans_type), 'axes': axes}
    for key, value in kwargs.items():
        transdict[unicode(key)] = value
    return transdict


def construct_reduce_dict(reduction_type, **kwargs):
    """Helper function used to construct reduction dictionary.

       Required Keyword arguments:

       :param reduction_type: The reduction type to be used.

       Optional Keyword arguments:

       :param axis: Some reduction functions support the axis parameters,
                    describing along which axis the reduction should be performed.
       :param x1: By default the reductions are performed on the output of the
                  previous data operation (x1='data'). We may reference the output
                  of, e.g., the fifth data operation by setting x1='data5'. x1
                  itself may also specify a separate data transformation and
                  reduction pipeline that operates on 'data'.
       :param min_dim: Minimum number of dimensions the input data should have in
                  order for the reduction should be applied.

       :returns: Dictionary with the description of the reduction operation.

    """
    if reduction_type not in reduction_allowed_functions:
        raise ValueError(unicode("Requested reduction operation not supported. Allowed operations are") +
                         unicode(reduction_allowed_functions.keys()))
    reducdict = {'reduction': unicode(reduction_type)}
    for key, value in kwargs.items():
        reducdict[unicode(key)] = value
    return reducdict


def transform_reduce_description_to_json(*args):
    """Convert the dictionary describing the transformation/reduction operations
       to a JSON string.

       :param args: The list or dictionaries with the description of the transformation
                     and reduction operations.

       :returns: JSON string
    """
    try:
        import json
    except ImportError:
        from django.utils import simplejson as json
    if len(args) == 1 and isinstance(args[0], list):
        return json.dumps(args[0])
    else:
        return json.dumps(construct_transform_reduce_list(*args))


def json_to_transform_reduce_description(json_string):
    """Convert the json string to the transformation/reduction dict.

       :param json_string: The json string to be converted.

       :returns: Python list or dict with the description
    """
    try:
        import json
    except ImportError:
        from django.utils import simplejson as json
    checkstring = is_transform_or_reduce(json_string)
    if not checkstring:
        warnings.warn(
            "The given JSON string may not define a valid transformation/reduction description.")
    transreduce = json.loads(json_string)
    if isinstance(transreduce, dict):
        return [transreduce]
    elif isinstance(transreduce, list):
        return transreduce
    else:
        if checkstring:
            warnings.warn(
                'The given JSON string may not define a valid transformation/reduction description')
        return transreduce

# def plot_transform_reduce_description(description):
#
#     import networkx as nx
#     import matplotlib.pyplot as plt
#
#     if isinstance(description,str) or isinstance(description,unicode):
#         alloperations = json_to_transform_reduce_description(description)
#     else:
#         alloperations = description
#
#     G = nx.Graph()
#     G.add_node('data')
#
#     def add_graph_node( operation, graph ):
#         if isinstance(operation, list):
#             nodenames = []
#             for co in operation:
#                 nodenames.append( add_graph_node(co, graph) )
#
#         elif isinstance(operation,dict):
#             try:
#                 ci = graph.add_node(operation['transformation'] , operation)
#             except:
#                 ci = graph.add_node(operation['reduction'] , operation)
#             for key, value in operation.items():
#                 if is_transform_or_reduce(value):
#                     cj = add_graph_node(value, graph)
#                     graph.add_edge(ci,cj)
#             return ci
#
#
#     add_graph_node( alloperations , G )
#     graph_pos = nx.spring_layout(G)
#     nx.draw(G, graph_pos)
#     plt.show()
