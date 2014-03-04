from omsi.shared.omsi_data_selection import *
import numpy as np

def main(argv=None):
    """Then main function"""

    npV1 = normalize_by_percentiles_V1()
    print ""
    print ""

    npV2 = normalize_by_percentiles_V2()
    print ""
    print ""

    print "#################################################################"
    print "##   Generating test data                                      ##"
    print "#################################################################"
    data0 = np.arange(25).reshape(5,5)
    print "Test data:"
    print data0
    print "#################################################################"
    print "##   Testing Normalize by percentile; Version 1                ##"
    print "#################################################################"
    npV1_secondary = {}
    npV1_result = transform_and_reduce_data(data=data0, operations=npV1, secondary_data=npV1_secondary)
    print "Result: Normalize by percentile; Version 1:"
    print npV1_result
    print ""
    print "#################################################################"
    print "##   Testing Normalize by percentile; Version 2                ##"
    print "#################################################################"
    npV2_secondary = {}
    npV2_result = transform_and_reduce_data(data=data0, operations=npV2, secondary_data=npV2_secondary)
    if np.all(npV1_result == npV2_result) :
        print "Results of normalization version 1 and 2 match!"
    else :
        print "Result: Normalize by percentile; Version 2:"
        print npV2_result
    #printDict(npV2_secondary)

def printDict(d):
    """Helper function used to print dicts in a more readable way"""
    print "____________________________________________"
    for key,value in d.items():
        print key
        print value
    print "____________________________________________"



def normalize_by_percentiles_V1( lower_percentile=5. , upper_percentile=95. ):
    """Define JSON string for normalization by: i) substracting a lower percentile,
       ii) dividing by an upper percentile, and iii) clipping the data so that
       values are in the range of 0 to 1.

       :param lower_percentile: The lower percentile to substract from the data (default=5)
       :param upper_percentile: The upper percentile the data should be divided by (default=95)

       :returns: JSON string with the description of the procedure
    """

    #1.1) Compute the lower percentile
    plow = construct_reduce_dict( reduction_type='percentile' ,
                                  q=lower_percentile )
    #1.2) Substract the lower percentile
    minus_p5 = construct_transform_dict( trans_type='arithmetic' ,
                                         operation='subtract',
                                         x1='data',
                                         x2=plow )
    #2.1) Compute the upper percentile
    phigh = construct_reduce_dict( reduction_type='percentile',
                                   q=upper_percentile )
    #2.2) Convert phigh to float
    phigh_as_float = construct_transform_dict( trans_type='astype' , dtype='float', x1=phigh, variable='phigh_as_float')
    #2.2) Divide by the upper percentile
    divide_by_p95 = construct_transform_dict( trans_type='arithmetic',
                                              operation='divide',
                                              x2=phigh_as_float )
    #3) Clip the data so that all values >1 are set to 1 and all values smaller than 0 are set to 0
    clip = construct_transform_dict( trans_type='scale',
                                     operation='clip',
                                     a_min=0,
                                     a_max=1 )
    operations_json = transform_reduce_description_to_json( minus_p5, divide_by_p95 , clip )

    print "#################################################################"
    print "##   Normalize by percentiles: Version 1                       ##"
    print "#################################################################"
    print "1) Purpose:"
    print ""
    print "Normalize the data by: i) subtracting by the "+str(lower_percentile)+"'tile,"
    print "                       ii) dividing by the "+str(upper_percentile)+"%'tile,"
    print "                       iii) clipping the data so that all values are in the range of [0,1]."
    print ""
    print "2) Illustration of the calculation performed:"
    print ""
    print "     /-->percentile(5)--\        /-->percentile(95)-->astype('float')--\ "
    print "    /                    \      /                                       \ "
    print "data------------->minus( -- )-->------------------------------->divide( -- )-->clip(0,1)-->output"
    print ""
    print "3) JSON description of the above calculation:"
    print ""
    print operations_json

    return operations_json


def normalize_by_percentiles_V2( lower_percentile=5. , upper_percentile=95. ):
    """Define JSON string for normalization by: i) substracting a lower percentile,
       ii) dividing by an upper percentile, and iii) clipping the data so that
       values are in the range of 0 to 1.

       :param lower_percentile: The lower percentile to substract from the data (default=5)
       :param upper_percentile: The upper percentile the data should be divided by (default=95)

       :returns: JSON string with the description of the procedure
    """

    #1) Compute the lower percentile from the input data
    plow = construct_reduce_dict( reduction_type='percentile' ,
                                  q=lower_percentile ,
                                  x1='data',
                                  variable='plow')
    #2) Compute the upper percentile from the input data
    phigh = construct_reduce_dict( reduction_type='percentile',
                                   q=upper_percentile,
                                   x1='data')
    phigh_minus_plow = construct_transform_dict( trans_type='arithmetic',
                                                 operation='subtract',
                                                 x2 = 'plow')
    phigh_as_float = construct_transform_dict( trans_type='astype',
                                               dtype='float',
                                               variable='phigh_as_float')

    #3) Substract the lower percentile --> Generates data3
    minus_plow = construct_transform_dict( trans_type='arithmetic' ,
                                         operation='subtract',
                                         x1='data',
                                         x2='plow')

    #4) Divide by the upper percentile --> Generates data4
    divide_by_phigh = construct_transform_dict( trans_type='arithmetic',
                                              operation='divide',
                                              x2='phigh_as_float')

    #3) Clip the data so that all values >1 are set to 1 and all values smaller than 0 are set to 0
    clip = construct_transform_dict( trans_type='scale',
                                     operation='clip',
                                     a_min=0,
                                     a_max=1)
    operations_json = transform_reduce_description_to_json( plow, phigh,
                                                            phigh_minus_plow,
                                                            phigh_as_float,
                                                            minus_plow,
                                                            divide_by_phigh ,
                                                            clip )

    print "######################################################################"
    print "##   Normalize by percentiles: Version 2 (explicit data references) ##"
    print "######################################################################"
    print "1) Purpose:"
    print ""
    print "Normalize the data by: i) subtracting by the "+str(lower_percentile)+"'tile,"
    print "                       ii) dividing by the "+str(upper_percentile)+"%'tile from the data before subtraction,"
    print "                       iii) clipping the data so that all values are in the range of [0,1]."
    print ""
    print "2) Illustration of the calculation performed:"
    print ""
    print "         /-->percentile(data,95)-->minus(..,plow)-->astype(..,'float')-->phigh_as_float-->\                                          \ "
    print "        /                                    /                                             \ "
    print "       /                             /------/                                               \ "
    print "      /-->percentile(data,5)-->plow--\                                                       \ "
    print "     /                                \                                                       \ "
    print "data--------------------->minus(data,plow)------------------------------->divide( .. , phigh_as_float )-->clip(0,1)-->output"
    print ""
    print "3) JSON description of the above calculation:"
    print ""
    print operations_json

    return operations_json




if __name__ == "__main__":
    main()