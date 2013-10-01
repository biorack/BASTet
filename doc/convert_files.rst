Converting and Accessing Files
==============================

.. code-block:: none
    
    ssh edison.nersc.gov
    cd /project/projectdirs/openmsi/devel/convert
    source setupEnvironment.csh
    python convertToOMSI.py <infile1> <output HDF5File>

Note if you use the bash shell then use ``setupEnvironment.bash`` instead. In order
to view a complete list of conversions options use: 

.. code-block:: bash

    python convertToOMSI.py --help

Making a converted file accesible to OpenMSI (Public)
-----------------------------------------------------

.. code-block:: bash

    ssh portal-auth.nersc.gov
    #<copy the converted file to /project/projectdirs/openmsi/omsi_data"
    cd /project/projectdirs/openmsi/omsi_data
    source update_filepermissions

Making a converted file accesilbe to OpenMSI (Private)
------------------------------------------------------

.. code-block:: bash

    ssh portal-auth.nersc.gov
    #copy the converted file to /project/projectdirs/openmsi/omsi_data_private"
    cd /project/projectdirs/openmsi/omsi_data_private
    source update_private_filepermissions

