# infile = "/nfs/disk1/users/amandabacon/LAPPD/lappd_analysis/output_lappd_0.root"
# infile = "/nfs/disk1/users/bharris/lappd_analysis/ppo-test_lappd_0.root"


import uproot
import awkward as ak
import matplotlib.pyplot as plt
import numpy as np

paths_scint = [
    "C:\\Users\\Evan\\Downloads\\lappd-158-lab1.7ppo-beta-no-mask-band-blocker-500nm-lp-centered-redo_lappd_1.root",
    "C:\\Users\\Evan\\Downloads\\lappd-158-lab1.7ppo-beta-mask-band-blocker-500nm-lp-centered_lappd_0.root",
    "C:\\Users\\Evan\\Downloads\\lappd-158-lab1.7ppo-beta-mask-band-blocker-500nm-lp-centered_lappd_1.root",
    "C:\\Users\\Evan\\Downloads\\lappd-158-lab1.7ppo-beta-no-mask-band-blocker-500nm-lp-centered-redo_lappd_0.root"
]

paths_water = [ "C:\\Users\\Evan\\Downloads\\lappd-158-water-beta-no-mask-band-blocker-500nm-lp-centered_lappd_1.root",
    "C:\\Users\\Evan\\Downloads\\lappd-158-water-beta-no-mask-band-blocker-500nm-lp-centered_lappd_0.root",
    "C:\\Users\\Evan\\Downloads\\lappd-158-water-beta-no-mask-500nm-lp-centered_lappd_0.root",
    "C:\\Users\\Evan\\Downloads\\lappd-158-water-beta-no-mask-500nm-lp-centered_lappd_1.root",
    "C:\\Users\\Evan\\Downloads\\lappd-158-water-beta-mask-500nm-lp-centered_lappd_0.root",
    "C:\\Users\\Evan\\Downloads\\lappd-158-water-beta-mask-500nm-lp-centered_lappd_1.root", ]



def max_calc(file_path):
    
    max = []
    n_tup = uproot.open(file_path)
    output = n_tup['output;1']
    event_info = output.arrays(output.keys(),library = "ak")
    waveforms = event_info['waveforms']
    lcns = event_info["lcns"]
    for i in range(len(waveforms)):
        for j in range(len(waveforms[i])):
            if(lcns[i][j] < 31):
                curr_max = ak.max(waveforms[i][j][10:-10])
                curr_ped = ak.mean(waveforms[i][j][500:900])
                curr_amp = curr_max - curr_ped
                max.append(curr_amp)
            
    return max

all_flat = max_calc(paths_scint[0])+ max_calc(paths_scint[1]) + max_calc(paths_scint[2]) + max_calc(paths_scint[3]) + max_calc(paths_water[0]) + max_calc(paths_water[1]) + max_calc(paths_water[2]) + max_calc(paths_water[3]) + max_calc(paths_water[4]) + max_calc(paths_water[5])
plt.hist(all_flat, bins = 100, range = (0, 300))


x_line = 55
plt.axvline(x=x_line, color="red", linestyle="--", linewidth=1)






plt.xlabel("amp in adc")
plt.title("Amplitudes of LAPPD Channels -- 55mV threshold")
plt.yscale("log")
plt.ylabel("counts")
plt.show()