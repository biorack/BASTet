 def memory_usage_resource(self):
     import sys
     import os
     import subprocess
     out = subprocess.Popen(['ps', 'v', '-p', str(os.getpid())],
     stdout=subprocess.PIPE).communicate()[0].split(b'\n')
     vsz_index = out[0].split().index(b'RSS')
     mem = float(out[1].split()[vsz_index]) / 1024
     return  mem
