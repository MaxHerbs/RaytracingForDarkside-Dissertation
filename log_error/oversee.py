from time import sleep
import os 


cmd = "scp godzilla:/storage/physics/phuddg/intensity_w_position/0.22/output_data/* output_data/"
while(True):
    complete = os.system(cmd)
    if complete == 0:
        print("Aquired Data")
        complete = os.system("python analysis.py")
    else:
        print("Failed to aquire data")
    sleep(180)