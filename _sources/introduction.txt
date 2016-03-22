Introduction
============

BASTet has been orginially developed as the analysis backend for the [OpenMSI](https://openmsi.nersc.gov/) science gateway. 

A central goal of BASTet is to facilitate shareable and reproducible analysis and to bridge the gaps between the various stages in the typical life-cycle of new analysis methods, e.g., when transition from a research prototype to production or when integrating an analysis into more complex workflows or when sharing analysis results. BASTet provides users with an environment that makes it easy to develop and deploy analyses via OpenMSI by providing: 

	* standardized analysis interfaces that make analyses easily accessible to users and enable developers to build an accessible ecosystem of domain analytics, 
	* automatic provenance for reproducible analytics,
	* standardized interfaces and a common format for storage, access, sharing and reuse of analysis results and raw MSI data 
	* support for integration of analyses to complex workflows, providing users a lightweight  and easy-to-use entry-point to define and manage analysis workflows, and
     * integrated tools and analytics to facilitate development and deployment. 


While the design of BASTet has been motivated by the needs of OpenMSI, its core design and functionality are much more broadly applicable to other applications as well.  Basket is implemented in Python using `NumPy` for data processing, `h5py` for HDF5-based data storage, and `mpi4py` for distributed parallel data processing. 


