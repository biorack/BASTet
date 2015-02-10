def mem_usage_system_psutil(self):
     """
     Memory usage as measured by the system's psutil tool.
     This can be usefule if the PYTHON psutil is not available.
     CAUTION: This function spawns a subprocess which can be
     problematic when having to measure memory usage many times.
     """
     import os
     import subprocess
     output = subprocess.Popen(['ps', 'v', '-p', str(os.getpid())],
                                stdout=subprocess.PIPE).communicate()[0].split(b'\n')
     vsz_index = output[0].split().index(b'RSS')
     mem = float(output[1].split()[vsz_index]) / 1024
     return  mem


def mem_usage_psutil():
    """
    Measure memory usage using python psutil.
    """
    import psutil
    p = psutil.Process(os.getpid())
    mem = p.get_memory_info()[0] / float(1024 * 1024)
    return mem


def mem_usage_resource():
    """
    Measure memory usage using python resource.
    While this approach can be faster then mem_usage_psutil(..)
    the resources module may not excluded orphand arrays deleted
    by the python interpreter are not always seen as cleared.
    """
    import resource
    rusage_denom = 1024.
    if sys.platform == 'darwin':
        # ... it seems that in OSX the output is different units ...
        rusage_denom = rusage_denom * rusage_denom
    mem = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / rusage_denom
    return mem
