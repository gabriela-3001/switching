import subprocess
import time
import os

cmd = 'bsub -q new-short -R "rusage[mem=500]" -u /dev/null ./a.out '

files = os.listdir()

for param in range(1,10202):
  if str(param) not in files:
    print(param)
    paramstr = ' '.join(map(str,[param]))
    _cmd = ' '.join([cmd, paramstr])
    subprocess.run(_cmd, shell=True)
    #print(cmd)
    time.sleep(0.005)