from omsi.analysis.findpeaks.experimental.omsi_lpf import omsi_lpf
from omsi.analysis.findpeaks.experimental.omsi_npg import omsi_npg
from omsi.analysis.findpeaks.experimental.omsi_peakcube import omsi_peakcube
from omsi.dataformat.omsi_file.main_file import omsi_file

import datetime
from time import time, ctime
import numpy as np
import os, sys, getopt, subprocess, tempfile
import string
from sys import argv,exit

def main(argv):

	totaltimekeeper = time()
	thisfilename = os.path.basename(__file__)

	# --- Default parameters ---

	# omsi file parameters
	omsiInFile = None
	expIndex = 0
	dataIndex = 0
	LPFIndex = None
	NPGIndex = None

	# skip analysis parameters
	SkipNPG = None
	SkipPeakCube = None

	# Flag for enabling peak cube through Carver big memory job
	QueuePeakCube = False

	# LPF default parameters
	myPeakheight = 10
	mySlwindow = 100
	mySmoothwidth = 3

	# NPG default parameters
	myMz_threshold = 0.05
	myTreecut = 0.1

	# repo parameter
	myRepo = None

	omsiparameters = [ "file=", "expIndex=", "dataIndex=","LPFIndex=", "NPGIndex=" ]
	skipparameters = [ "SkipNPG", "SkipPeakCube" ]
	lpfparameters = [ "peakheight=", "slwindow=", "smoothwidth=" ]
	npgparameters = [ "mz_threshold=", "treecut="]
	extraparameters = ["QueuePeakCube", "repo="]

	myparameters = omsiparameters + skipparameters + lpfparameters + npgparameters + extraparameters

	helpstring = "Call \"python " + thisfilename + " -h\" for help information."

	if(len(argv) < 1):
		print "Usage error! " + helpstring
		exit(0)

	try:
		opts, args = getopt.getopt(argv,"hf:", myparameters)
	except getopt.GetoptError:
		print "Parameter Error: " + helpstring
		sys.exit(2)

	for opt, arg in opts:
		if opt == '-h':
			printHelp(thisfilename)
			sys.exit()
		try:
			if opt in ("-f", "--file"):
				omsiInFile = arg
			elif opt == "--expIndex":
				expIndex = int(arg)
			elif opt == "--dataIndex":
				dataIndex = int(arg)
			elif opt == "--LPFIndex":
				LPFIndex = int(arg)
			elif opt == "--NPGIndex":
				NPGIndex = int(arg)
			elif opt == "--SkipNPG":
				SkipNPG = 1
			elif opt == "--SkipPeakCube":
				SkipPeakCube = 1
			elif opt == "--peakheight":
				myPeakheight = int(arg)
			elif opt == "--slwindow":
				mySlwindow = int(arg)
			elif opt == "--smoothwidth":
				mySmoothwidth = int(arg)
			elif opt == "--mz_threshold":
				myMz_threshold = float(arg)
			elif opt == "--treecut":
				myTreecut = float(arg)
			elif opt == "--QueuePeakCube":
				QueuePeakCube = True
			elif opt == "--repo":
				myRepo = arg
		except:
			print "Parameter Error: " + helpstring
			exit(0)

	# Open the input HDF5 file
	if omsiInFile is None:
		print "No OMSI file provided! " + helpstring
		exit(0)

	print "Input file: " , omsiInFile
	print "OMSI parameters: expIndex =", expIndex, ", dataIndex =", dataIndex
	print "Analysis indexes: LPFIndex =", LPFIndex, ", NPGIndex =", NPGIndex
	sys.stdout.flush()

	# run flags
	runlpf = 0
	runnpg = 0
	runpc = 0

	# --- LPF
	if(LPFIndex is None):
		runlpf = 1
		LPFanalysisindex = run_lpf(omsiInFile, expIndex, dataIndex, ph=myPeakheight, slw=mySlwindow, smw=mySmoothwidth)
	else:
		LPFanalysisindex = LPFIndex

	# --- NPG
	if(NPGIndex is None and SkipNPG is None):
		runnpg = 1
		NPGanalysisindex = run_npg(omsiInFile, expIndex, dataIndex, LPFanalysisindex, mzth=myMz_threshold, tcut=myTreecut)
	elif(NPGIndex is not None):
		NPGanalysisindex = NPGIndex
	else:
		NPGanalysisindex = None

	# run peak cube in big memory node
	if(QueuePeakCube == True and LPFanalysisindex is not None and NPGanalysisindex is not None):
		pcstring = "python " + thisfilename + " -f " + omsiInFile + " --expIndex=" + str(expIndex) + " --dataIndex=" + str(dataIndex)
		pcstring += " --LPFIndex=" + str(LPFanalysisindex) + " --NPGIndex=" + str(NPGanalysisindex) + " \n"

		jobname = queuePCjob(pcstring, therepo = myRepo)
		SkipPeakCube = 1

		print "Peak Cube job name:", jobname
		print "Peak Cube queued."
		print "Peak Cube output file: pc_job." + jobname + ".out.txt"

	# --- Peak Cube
	if(SkipPeakCube is None and ((NPGIndex is None and SkipNPG is None) or NPGIndex is not None)):
		runpc = 1
		PCanalysisindex = run_peakcube(omsiInFile, expIndex, dataIndex, LPFanalysisindex, NPGanalysisindex)
	else:
		PCanalysisindex = None

	print "--- all complete ---"
	print "LPF Analysis Index:", LPFanalysisindex
	print "NPG Analysis Index:", NPGanalysisindex
	print "PC  Analysis Index:", PCanalysisindex
	print "Input file: " , omsiInFile
	print "Total time: ", round((time() - totaltimekeeper) / 60, 2), "min"

# queue peak cube job
def queuePCjob(pcstring, therepo = None):

	# create temporary pbs script for carver
	(sfd, scriptfile) = tempfile.mkstemp(prefix="pcscript.", suffix=".pbs")
	pbsCreation = generateScript(scriptfile, PFcontent = pcstring, repo=therepo)
	os.close(sfd)

	if not pbsCreation:
		print "Error creating pbs script."
		exit(0)

	print "\nQueueing carver peak cube job..."

	qsubout = subprocess.check_output(["qsub", scriptfile])
	myjobid = string.split(qsubout, ".")[0]
	myjobname = string.split(qsubout, "\n")[0]

	print "Job id:", myjobid

	# delete temporary scriptfile
	os.remove(scriptfile)

	return myjobname

# lpf analysis
def run_lpf(omsiInFile, expIndex, dataIndex, ph, slw, smw):

	#Open the input HDF5 file
	try:
		omsiFile = omsi_file(omsiInFile,'r+')
	except:
		print "Error opening input file \"",omsiInFile,"\": ", sys.exc_info()[0]
		exit(0)


	# Get the experiment and data
	exp = omsiFile.get_experiment( expIndex )
	data = exp.get_msidata(dataIndex)
	peaksMZdata = data.mz[:]

	# LPF --------------
	print "\n--- Executing LPF ---"
	myLPF = omsi_lpf(name_key="omsi_lpf_" + str(ctime()))
	myLPF.execute( msidata=data, mzdata=peaksMZdata, peakheight = ph, slwindow = slw, smoothwidth = smw)
	print "\n\nResults"
	peaksBins = myLPF['LPF_Peaks_MZ'][:]
	print "peaksBins:\n", peaksBins
	peaksIntensities = myLPF['LPF_Peaks_Vals'][:]
	print "peaksIntensities:\n", peaksIntensities
	peaksArrayIndex = myLPF['LPF_Peaks_ArrayIndex'][:]
	print "peaksArrayIndex:\n", peaksArrayIndex

	print "\nSaving HDF5 analysis..."
	LPFanalysis, LPFanalysisindex = exp.create_analysis(myLPF)
	print "done!"
	print "--- omsi_lpf complete ---\n"
	print "LPF Analysis Index:", LPFanalysisindex

	# flush of lpf
	exp.experiment.file.flush()
	omsiFile.close_file()

	return LPFanalysisindex

# npg analysis
def run_npg(omsiInFile, expIndex, dataIndex, LPFIndex, mzth, tcut):

	#Open the input HDF5 file
	try:
		omsiFile = omsi_file(omsiInFile,'r+')
	except:
		print "Error opening input file \"",omsiInFile,"\": ", sys.exc_info()[0]
		exit(0)

	# Get the experiment and data
	exp = omsiFile.get_experiment( expIndex )
	data = exp.get_msidata(dataIndex)
	LPFanalysis = exp.get_analysis( LPFIndex )

	peaksBins = LPFanalysis['LPF_Peaks_MZ'][:]
	peaksIntensities = LPFanalysis['LPF_Peaks_Vals'][:]
	peaksArrayIndex = LPFanalysis['LPF_Peaks_ArrayIndex'][:]
	peaksMZdata = data.mz[:]
	peaksMZ = peaksMZdata[peaksBins]

	# NPG --------------
	print "\n--- Executing NPG ---"
	myNPG = omsi_npg(name_key="omsi_npg_" + str(ctime()))
	#myNPG.omsi_npg_exec(peaksBins, peaksIntensities, peaksArrayIndex, peaksMZdata, peaksMZ, MZ_TH = mzth, clusterCut = tcut)
	myNPG.execute( peaksBins=peaksBins,
                   npg_peaks_Intensities=peaksIntensities,
                   npg_peaks_ArrayIndex=peaksArrayIndex,
                   peaksMZdata=peaksMZdata,
                   peaksMZ=peaksMZ,
                   npg_mz_threshold=mzth,
                   npg_cluster_treecut=tcut   )
	print "\n\nResults"
	NPGPL = myNPG['npghc_peaks_labels']
	print "NPG HC Peaks Labels: \n", NPGPL
	NPGLL = myNPG['npghc_labels_list']
	print "NPG HC Labels List: \n", NPGLL

	print "\nSaving HDF5 analysis..."
	NPGanalysis, NPGanalysisindex = exp.create_analysis(myNPG)
	print "done!"
	print "--- omsi_npg complete ---\n"
	print "NPG Analysis Index:", NPGanalysisindex

	# flush of npg
	exp.experiment.file.flush()
	omsiFile.close_file()

	return NPGanalysisindex

# peak cube generation
def run_peakcube(omsiInFile, expIndex, dataIndex, LPFIndex, NPGIndex):

	#Open the input HDF5 file
	try:
		omsiFile = omsi_file(omsiInFile,'r+')
	except:
		print "Error opening input file \"",omsiInFile,"\": ", sys.exc_info()[0]
		exit(0)

	# Get the experiment
	exp = omsiFile.get_experiment( expIndex )
	data = exp.get_msidata(dataIndex)
	peaksMZdata = data.mz[:]
	LPFanalysis = exp.get_analysis( LPFIndex )
	NPGanalysis = exp.get_analysis( NPGIndex )

	peaksBins = LPFanalysis['LPF_Peaks_MZ'][:]
	peaksIntensities = LPFanalysis['LPF_Peaks_Vals'][:]
	peaksArrayIndex = LPFanalysis['LPF_Peaks_ArrayIndex'][:]
	NPGPeaksLabels = NPGanalysis['npghc_peaks_labels'][:]
	NPGLabelsList = NPGanalysis['npghc_labels_list'][:]

	print "\n--- Creating Peak Cube ---"
	myPC = omsi_peakcube(name_key="omsi_peakcube_" + str(ctime()))
	#myPC.omsi_peakcube_exec(peaksBins, peaksIntensities, peaksArrayIndex, peaksMZdata, NPGPeaksLabels, NPGLabelsList)
	myPC.execute( peaksBins = peaksBins,
                  peaksIntensities = peaksIntensities,
                  peaksArrayIndex = peaksArrayIndex,
                  peaksMZdata = peaksMZdata,
                  HCpeaksLabels = NPGPeaksLabels ,
                  HCLabelsList = NPGLabelsList)


	PCm = myPC['npg_peak_cube_mz']
	print "NPG Peak Cube Mzs: \n", PCm.shape, "\n", PCm
	PMz = myPC['npg_peak_mz']
	print "NPG Peak Mz: \n", PMz.shape, "\n", PMz

	print "\nSaving HDF5 analysis..."
	PCanalysis, PCanalysisindex = exp.create_analysis(myPC)
	print "done!"
	print "--- omsi_peakcube complete ---\n"
	print "PeakCube Analysis Index:", PCanalysisindex


	# flush of peakcube
	exp.experiment.file.flush()
	omsiFile.close_file()

	return PCanalysisindex

# generates pfscript.pbs for carver big memory job
def generateScript(scriptfile, PFcontent = None, repo = None):

	mypbs  = "#PBS -q serial\n"
	mypbs += "#PBS -l pvmem=40GB \n"
	mypbs += "#PBS -l walltime=06:00:00 \n"
	if(repo is not None):
		mypbs += "#PBS -A " + repo +" \n"
	mypbs += "#PBS -N pc_job \n"
	mypbs += "#PBS -e pc_job.$PBS_JOBID.err.txt \n"
	mypbs += "#PBS -o pc_job.$PBS_JOBID.out.txt \n"
	mypbs += "source /global/project/projectdirs/openmsi/sources/peakfinding_env.txt \n"
	mypbs += "cd $PBS_O_WORKDIR \n"
	mypbs += "echo '[peak cube start]' \n"
	mypbs += "date \n"

	if(PFcontent):
		mypbs += "echo '+ pc run: ' \n"
		mypbs += PFcontent
		mypbs += "echo '+ pc complete' \n"
		mypbs += "date \n"

	mypbs += "echo '[peak cube end]' \n"
	mypbs += "date \n"

	if(PFcontent is None):
		scriptcheck = 0
	else:
		scriptcheck = 1

	try:
		pbs_file = open(scriptfile, "w+")
		pbs_file.write(mypbs)
		pbs_file.close()
	except:
		scriptcheck = 0

	return scriptcheck

# help documentation
def printHelp(thisfilename):
	print "\n---------- Peak Finder Help: ----------"
	print "To run the analysis use \"-f OMSIFile\" or \"--file=OMSIFile\", for example:"
	print "\tpython " + thisfilename + " -f OMSIFile"
	print "\nThis will run local and global peak finding on a particular OMSIFile"
	print "with the default omsi file parameters:"
	print "\texpIndex=0, dataIndex=0"
	print "and the default local peak finding parameters:"
	print "\tpeakheight=10, slwindow=100, smoothwidth=3"
	print "and the default global peak finding parameters:"
	print "\tmz_threshold=0.05, treecut=0.1"
	print "\nIf different parameter values are desired please specify as arguments"
	print "preceded by \"--\", for example: "
	print "\n\tpython " + thisfilename + " -f OMSIFile --peakheight=15 --smoothwidth=5"
	print "\nNote: It is recommended to provide an appropriate peakheight parameter"
	print "corresponding to the data properties."
	print "\nThe program runs both LPF and NPG analysis and PeakCube generation by default."
	print "Note that NPG analysis depends on LPF results and PeakCube generation depends"
	print "on NPG results. If desired, one can run a specific analysis by providing the"
	print "analysis index of its dependency in the parameters, with \"LPFIndex\" and \"NPGIndex\"."
	print "For example, to run just the NPG analysis and PeakCube generation you can provide"
	print "the index of the LPF analysis and call the program with:"
	print "\n\tpython " + thisfilename + " -f OMSIFile --LPFIndex=1"
	print "\nThis will run the NPG analysis and PeakCube generation, skipping the LPF"
	print "analysis and reading the necessary LPF data from analysis index 1 in OMSIFile."
	print "Similarly, to run only the PeakCube generation one can call the program with:"
	print "\n\tpython " + thisfilename + " -f OMSIFile --LPFIndex=1 --NPGIndex=2"
	print "\nTo skip a particular section of the analysis one can use the parameters \"SkipNPG\""
	print "and \"SkipPeakCube\". For example, to run just the LPF analysis without NPG,"
	print "one can call:"
	print "\n\tpython " + thisfilename + " -f OMSIFile --SkipNPG"
	print "\nOr to run both LPF and NPG but not the PeakCube generation, call:"
	print "\n\tpython " + thisfilename + " -f OMSIFile --SkipPeakCube"
	print "\nNote that PeakCube generation not be performed if \"SkipNPG\" is enabled. If this"
	print "is not desired then use the \"NPGIndex\" parameter instead."
	print "\nBy default, the PeakCube generation is done interactively. Optionally, to queue"
	print "it as a carver big memory job use the parameter 'QueuePeakCube'. For example:"
	print "\n\tpython " + thisfilename + " -f OMSIFile --peakheight=15 --QueuePeakCube"
	print "\nThis will run the LPF and NPG analysis and then queue the PeakCube generation."
	print "To see the output of the job open the indicated file with prefix 'pc_job'. For"
	print "example:"
	print "\n\t'pc_job.JOBID.cvrsvc09-ib.out.txt'"
	print "\nThe output file will be available once the job is complete."
	print "---------------------------------------"

# stop python
def stop():
	raw_input("Stop!")

if __name__ == "__main__":
	main(sys.argv[1:])
