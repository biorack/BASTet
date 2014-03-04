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
    print "##   Testing the different functions                           ##"
    print "#################################################################"
    data0 = np.arange(25).reshape(5,5)
    print "Test data:"
    print data0
    npV1_result = transform_and_reduce_data(data=data0, operations=npV1)
    print npV1_result
    print "Result: Normalize by percentile; Version 1:"
    print "Test data:"
    print data0
    npV2_result = transform_and_reduce_data(data=data0, operations=npV2)
    print "Result: Normalize by percentile; Version 2:"
    print npV2_result




def normalize_by_percentiles_V1( lower_percentile=5 , upper_percentile=95 ):
    """Define JSON string for normalization by: i) substracting a lower percentile,
       ii) dividing by an upper percentile, and iii) clipping the data so that
       values are in the range of 0 to 1.

       :param lower_percentile: The lower percentile to substract from the data (default=5)
       :param upper_percentile: The upper percentile the data should be divided by (default=95)

       :returns: JSON string with the description of the procedure
    """

    #1.1) Compute the lower percentile
    plow = construct_reduce_dict( reduction_type='percentile' ,
                                  p=lower_percentile )
    #1.2) Substract the lower percentile
    minus_p5 = construct_transform_dict( trans_type='arithmetic' ,
                                         operation='subtract',
                                         x1='data',
                                         x2=plow )
    #2.1) Compute the upper percentile
    phigh = construct_reduce_dict( reduction_type='percentile',
                                 p=upper_percentile )
    #2.2) Convert data to float
    as_float = construct_transform_dict( trans_type='astype' , dtype='float')
    #2.3) Convert the percentile value to float to make sure that we get float data rather than integers
    phigh_as_float = construct_transform_reduce_list( phigh, as_float )
    #2.2) Divide by the upper percentile
    divide_by_p95 = construct_transform_dict( trans_type='arithmetic',
                                              operation='divide',
                                              x1='data',
                                              x2=phigh_as_float )
    #3) Clip the data so that all values >1 are set to 1 and all values smaller than 0 are set to 0
    clip = construct_transform_dict( trans_type='scale',
                                     operation='clip',
                                     a_min=0,
                                     a_max=1)
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
    print "      /-->percentile(5)--\             /-->percentile(95)-->astype('float')--\ "
    print "     /                    \           /                                       \ "
    print "data0------------->minus( -- )-->data1------------------------------->divide( -- )-->data2-->clip(0,1)-->output"
    print ""
    print "3) JSON description of the above calculation:"
    print ""
    print operations_json

    return operations_json


def normalize_by_percentiles_V2( lower_percentile=5 , upper_percentile=95 ):
    """Define JSON string for normalization by: i) substracting a lower percentile,
       ii) dividing by an upper percentile, and iii) clipping the data so that
       values are in the range of 0 to 1.

       :param lower_percentile: The lower percentile to substract from the data (default=5)
       :param upper_percentile: The upper percentile the data should be divided by (default=95)

       :returns: JSON string with the description of the procedure
    """

    #1) Compute the lower percentile from data0 --> Generates data1
    plow = construct_reduce_dict( reduction_type='percentile' ,
                                  p=lower_percentile ,
                                  x1='data0')
    #2) Compute the upper percentile from data0 and convert to float --> Generates data2
    phigh = construct_reduce_dict( reduction_type='percentile',
                                   p=upper_percentile,
                                   x1='data0')
    as_float = construct_transform_dict( trans_type='astype' , dtype='float')
    phigh_as_float = construct_transform_reduce_list( phigh, as_float )

    #3) Substract the lower percentile --> Generates data3
    minus_p5 = construct_transform_dict( trans_type='arithmetic' ,
                                         operation='subtract',
                                         x1='data0',
                                         x2='data1' )

    #4) Divide by the upper percentile --> Generates data4
    divide_by_p95 = construct_transform_dict( trans_type='arithmetic',
                                              operation='divide',
                                              x1='data3',
                                              x2='data2' )

    #3) Clip the data so that all values >1 are set to 1 and all values smaller than 0 are set to 0
    clip = construct_transform_dict( trans_type='scale',
                                     operation='clip',
                                     x1='data4',
                                     a_min=0,
                                     a_max=1)
    operations_json = transform_reduce_description_to_json( plow, phigh_as_float, minus_p5, divide_by_p95 , clip )

    print "######################################################################"
    print "##   Normalize by percentiles: Version 2 (explicit data references) ##"
    print "######################################################################"
    print "1) Purpose:"
    print ""
    print "Normalize the data by: i) subtracting by the "+str(lower_percentile)+"'tile,"
    print "                       ii) dividing by the "+str(upper_percentile)+"%'tile,"
    print "                       iii) clipping the data so that all values are in the range of [0,1]."
    print ""
    print "2) Illustration of the calculation performed:"
    print ""
    print "        /-->percentile(5)-->astype('float')-->data2-----\ "
    print "       /                                                 \ "
    print "      /-->percentile(5)-->data1--\                        \ "
    print "     /                            \                        \ "
    print "data0--------------------->minus( -- )-->data3---->divide( -- )-->data4-->clip(0,1)-->output"
    print ""
    print "3) JSON description of the above calculation:"
    print ""
    print operations_json

    return operations_json




if __name__ == "__main__":
    main()