import datetime
from time import time, ctime, sleep
import os, sys, getopt, subprocess, tempfile
import string
from sys import argv, exit


def main(argv):
    totaltimekeeper = time()

    thisfilename = os.path.basename(__file__)
    driverfilename = "pfrun.py"
    PFjobname = None
    PFjobid = None
    PCjobname = None
    PCjobid = None

    # omsi file parameters
    omsiInFile = None
    expIndex = 0
    dataIndex = 0
    LPFIndex = None
    NPGIndex = None

    # extra analysis parameters
    SkipNPG = None
    SkipPeakCube = None
    monitor = None

    # LPF default parameters
    myPeakheight = 10
    mySlwindow = 100
    mySmoothwidth = 3

    # NPG default parameters
    myMz_threshold = 0.05
    myTreecut = 0.1


    # repo parameter
    myRepo = None

    omsiparameters = ["file=", "expIndex=", "dataIndex=", "LPFIndex=", "NPGIndex="]
    extraparameters = ["SkipNPG", "SkipPeakCube", "monitor", "repo="]
    lpfparameters = ["peakheight=", "slwindow=", "smoothwidth="]
    npgparameters = ["mz_threshold=", "treecut="]

    myparameters = omsiparameters + extraparameters + lpfparameters + npgparameters

    helpstring = "Call \"python " + thisfilename + " -h\" for help information."

    if (len(argv) < 1):
        print "Usage error! " + helpstring
        exit(0)

    try:
        opts, args = getopt.getopt(argv, "hf:", myparameters)
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
            elif opt == "--monitor":
                monitor = 1
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
            elif opt == "--repo":
                myRepo = arg
        except:
            print "Parameter Error: " + helpstring
            exit(0)

    # Open the input HDF5 file
    if omsiInFile is None:
        print "No OMSI file provided! " + helpstring
        exit(0)
    else:
        print "OMSI file:", omsiInFile

    # python call strings
    pfstring = "python " + driverfilename + " -f " + omsiInFile + " --expIndex=" + str(
        expIndex) + " --dataIndex=" + str(dataIndex)

    if (myRepo is not None):
        pfstring += " --repo=" + myRepo

    lpfstring = " --peakheight=" + str(myPeakheight) + " --slwindow=" + str(mySlwindow) + " --smoothwidth=" + str(
        mySmoothwidth)
    npgstring = " --mz_threshold=" + str(myMz_threshold) + " --treecut=" + str(myTreecut)

    if (SkipPeakCube is None and ((NPGIndex is None and SkipNPG is None) or NPGIndex is not None)):
        runpc = True
    else:
        runpc = False

    # all in one
    PFcmd = getPFcmd(pfstring, lpfstring, npgstring, LPFIndex, NPGIndex, SkipNPG, SkipPeakCube)
    if not PFcmd:
        print "Parameter Error: " + helpstring
        exit(0)


    # create temporary pbs script
    (sfd, scriptfile) = tempfile.mkstemp(prefix="pfscript.", suffix=".pbs")
    pbsCreation = generateScript(scriptfile, PFcontent=PFcmd, repo=myRepo)
    os.close(sfd)

    if not pbsCreation:
        print "Error creating pbs script."
        exit(0)

    print "PBS script created:", scriptfile
    print "Queueing job..."

    qsubout = subprocess.check_output(["qsub", scriptfile])
    PFjobid = string.split(qsubout, ".")[0]
    PFjobname = string.split(qsubout, "\n")[0]

    print "Job id:", PFjobid

    # delete temporary scriptfile
    os.remove(scriptfile)

    if (monitor):
        print "\n[PF Monitoring start]"

        PCjobname = monitorJob(PFjobid, PFjobname, runtype="pf")

        print "\n[PF Monitoring complete]"

        if (PCjobname is not None):
            print "\n[Peak Cube Monitoring start]"
            PCjobid = string.split(PCjobname, ".")[0]
            PCresult = monitorJob(PCjobid, PCjobname, runtype="pc")
            print "\n[Peak Cube Monitoring complete]"

    if (PFjobname is not None or PCjobname is not None):
        if (PFjobname):
            print "Peak finding output file: pf_job." + PFjobname + ".out.txt"
        if (PCjobname):
            print "Peak cube output file: pc_job." + PCjobname + ".out.txt"
        elif (runpc):
            print "Peak cube output file: [To be queued] check 'pf_job." + PFjobname + ".out.txt' after peak finding is complete."

    print "--- fin ---"


# +++++++++++++++++++++++ helper functions ++++++++++++++++++++++++++
# monitors a job output
def monitorJob(jobid, jobname, runtype="pf"):
    runkeystring = "State: Running"
    completekeystring = "State: Completed"
    exitkeystring = "Jobs exit status"
    pcskeystring = "Peak Cube job name: "
    jobisrunning = False
    keepmonitoring = True
    currentline = 0
    checkjobsleep = 10
    filereadsleep = 5
    PCjobname = None
    ioerrorcount = 0
    ioerrormax = 3  # tolerance amount for error in reading output file. More attempts than this will stop monitoring

    print "Waiting for job to start...",
    sys.stdout.flush()

    while (jobisrunning == False):
        try:
            processoutput = subprocess.check_output(["checkjob", str(jobid)], stderr=subprocess.STDOUT)
            runresult1 = string.find(processoutput, runkeystring)
            runresult2 = string.find(processoutput, completekeystring)
        except:
            runresult1 = -1
            runresult2 = -1

        if (runresult1 != -1 or runresult2 != -1):
            jobisrunning = True
            print "\nJob is running..."
        else:
            sleep(checkjobsleep)

    print "--- [job output start] ---"

    while (keepmonitoring):

        joboutput = getJobOutput(jobname, runtype)

        if (joboutput is None):
            if (ioerrorcount > ioerrormax):
                print "Error reading outputfile for job:", jobid
                keepmonitoring = False
            else:
                ioerrorcount += 1
                sleep(2)
        else:
            lineamt = len(joboutput) - currentline

            for i in xrange(lineamt):
                thisline = joboutput[currentline]

                # check for incomplete analysis line
                if (i == lineamt - 1):
                    startlpfline = string.find(thisline, "smooth - done")
                    endlpfline = string.find(thisline, "peakdet - done")
                    if (startlpfline != -1 and endlpfline != -1):
                        print thisline,
                    elif (thisline[-1] == '\r'):
                        print thisline[0:-1] + '\r',
                        currentline -= 1
                    else:
                        currentline -= 1
                else:
                    print thisline,

                foundpc = string.find(thisline, pcskeystring)
                if (foundpc != -1):
                    PCjobname = string.split(thisline, pcskeystring)[1].strip()

                foundend = string.find(thisline, exitkeystring)
                if (foundend != -1):
                    keepmonitoring = False

                currentline += 1
                sys.stdout.flush()

            if (keepmonitoring):
                sleep(filereadsleep)

    print "\n--- [job output complete] ---"
    return PCjobname


# gets the job output
def getJobOutput(jobname, runtype):
    jobcomplete = False
    myoutput = None

    oufile = jobname + ".OU"
    if (runtype == "pf"):
        pffile = "pf_job." + jobname + ".out.txt"
    else:
        pffile = "pc_job." + jobname + ".out.txt"

    try:
        myfile = open(oufile, 'r')
        myoutput = myfile.readlines()
        myfile.close()
    except:
        myoutput = None
        jobcomplete = True

    if (jobcomplete):
        # in case of transitioning from running to complete,
        # check filename for a complete job after 2 seconds
        sleep(2)
        try:
            myfile = open(pffile, 'r')
            myoutput = myfile.readlines()
            myfile.close()
        except:
            myoutput = None

    return myoutput


# generates pfscript.pbs
def getPFcmd(pfstring, lpfstring, npgstring, LPFIndex, NPGIndex, SkipNPG, SkipPeakCube, IndexFile=None):
    runlpf = False
    runnpg = False
    runpc = False

    # set run flags
    if (LPFIndex is None):
        runlpf = True
    if (NPGIndex is None and SkipNPG is None):
        runnpg = True
    if (SkipPeakCube is None and ((NPGIndex is None and SkipNPG is None) or NPGIndex is not None)):
        runpc = True

    # builf pf cmd
    if (runlpf and runnpg and runpc):
        pfstring += lpfstring
        pfstring += npgstring
        PFcmd = pfstring
    elif (runlpf and runnpg):
        pfstring += lpfstring
        pfstring += npgstring
        PFcmd = pfstring + " --SkipPeakCube"
    elif (runlpf and runpc):
        pfstring += lpfstring
        PFcmd = pfstring + " --NPGIndex=" + str(NPGIndex)
    elif (runnpg and runpc):
        pfstring += npgstring
        PFcmd = pfstring + " --LPFIndex=" + str(LPFIndex)
    elif (runlpf):
        pfstring += lpfstring
        PFcmd = pfstring + " --SkipNPG --SkipPeakCube"
    elif (runnpg):
        pfstring += npgstring
        PFcmd = pfstring + " --LPFIndex=" + str(LPFIndex) + " --SkipPeakCube"
    elif (runpc):
        pfstring += " --LPFIndex=" + str(LPFIndex) + " --NPGIndex=" + str(NPGIndex)
        PFcmd = pfstring
    else:
        PFcmd = None

    if (PFcmd is not None):
        if (runpc):
            PFcmd += " --QueuePeakCube"

        PFcmd += " \n"

    return PFcmd


# generates pfscript.pbs for dirac
def generateScript(scriptfile, PFcontent=None, repo=None):
    mypbs = "#PBS -q dirac_reg\n"
    mypbs += "#PBS -l nodes=1:ppn=8:fermi \n"
    mypbs += "#PBS -l walltime=06:00:00 \n"
    if (repo is not None):
        mypbs += "#PBS -A " + repo + " \n"
    mypbs += "#PBS -N pf_job \n"
    mypbs += "#PBS -e pf_job.$PBS_JOBID.err.txt \n"
    mypbs += "#PBS -o pf_job.$PBS_JOBID.out.txt \n"
    mypbs += "echo '[peakfinder start]' \n"
    mypbs += "source /global/project/projectdirs/openmsi/sources/peakfinding_env.txt \n"
    mypbs += "cd $PBS_O_WORKDIR \n"
    mypbs += "date \n"

    if (PFcontent):
        mypbs += "echo '+ pf run: ' \n"
        mypbs += PFcontent
        mypbs += "echo '+ pf complete' \n"
        mypbs += "date \n"

    mypbs += "echo '[peakfinder end]' \n"
    mypbs += "date \n"

    if (PFcontent is None):
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
    print "\nThe output can be monitored by using the --monitor flag. This will read the"
    print "output files generated by dirac and carver jobs. For example:"
    print "\n\tpython " + thisfilename + " -f OMSIFile --peakheight=15 --monitor"
    print "\nIf the monitor flag is not used, the output will be written in the current working"
    print "directory with the prefix 'pf_job' for LPF and NPG analysis and 'pc_job' for peak"
    print "cube generation. For example:"
    print "\n\t'pf_job.JOBID.cvrsvc09-ib.out.txt' and 'pc_job.JOBID.cvrsvc09-ib.out.txt'"
    print "\nThis file queues a job in dirac for lpf and npg analysis and the dirac job can queue"
    print "a big memory job in carver if the peak cube will be generated. The monitoring can be"
    print "performed by reading each job output file. Temporary files are created to generate"
    print "customized pbs scripts and then removed once the job is queued."
    print "---------------------------------------"


# stop python
def stop():
    raw_input("Stop!")


if __name__ == "__main__":
    main(sys.argv[1:])
