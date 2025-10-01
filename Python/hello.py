import uproot
import random
import awkward as ak 
import numpy as np
import matplotlib.pyplot as plt


file_path1 = "C:\\Users\\Evan\\Downloads\\lappd-158-water-beta-mask-500nm-lp-centered_lappd_1.root"
file_path0 = "C:\\Users\\Evan\\Downloads\\lappd-158-water-beta-mask-500nm-lp-centered_lappd_0.root"

n_tup = uproot.open(file_path0)
output = n_tup['output;1']
event_info = output.arrays(output.keys(),library = "ak")
waveforms = event_info['waveforms']
peaks =event_info['peaks'] 
ped = event_info['pedestals'] 
lcns =  event_info['lcns']
n_abv_thresh = event_info['nsamples_above_threshold']
times = event_info['times']

n= 1 ## n above thresh mask 
amp_mask = 10 ## in Volts

def residualscalc(file_path):
    n_tup = uproot.open(file_path)
    output = n_tup['output;1']
    event_info = output.arrays(output.keys(),library = "ak")
    waveforms = event_info['waveforms']
    peaks =event_info['peaks'] 
    ped = event_info['pedestals'] 
    lcns =  event_info['lcns']
    n_abv_thresh = event_info['nsamples_above_threshold']
    times = event_info['times']
    t_res = []
    discard = []
    

    
    def locate_groputrigger_time(i, j):
        if lcns[i][j] < 8:
            pedestal = ped[i][j]
            wf = waveforms[i][31][100:-100]
            crossed = wf - pedestal > (0.6 * (peaks[i][j] - pedestal))
            above = ak.where(crossed)[0]
            return 0.2 * above[0]
        
        elif 7 < lcns[i][j] < 16:
            pedestal = ped[i][j]
            wf = waveforms[i][31][100:-100]
            crossed = wf - pedestal > (0.6 * (peaks[i][j] - pedestal))
            above = ak.where(crossed)[0]
            return 0.2 * above[0]
        
        elif 15 < lcns[i][j] < 24:
            return times[i][33]
        elif 23 < lcns[i][j] < 31:
            return times[i][34]
        
        
       
    for i in range(len(waveforms)):
        discard_channel_per_event = []
        t_res_event = []
        for j in range(len(waveforms[i])):
            wf = waveforms[i][j]
            amp = peaks[i][j] - ped[i][j]
            if (n_abv_thresh[i][j] > n) and (amp > amp_mask) :
                if lcns[i][j] not in [31,32,48,64,80]:
                    residual = times[i][j] - locate_groputrigger_time(i,j)
                    t_res_event.append(residual)
                else:
                     discard_channel_per_event.append(j)
            else : 
                discard_channel_per_event.append(j)

        discard.append(discard_channel_per_event)
        t_res.append(t_res_event) 

    t_res_arr = ak.Array(t_res)
    discard_arr = ak.Array(discard)
    return t_res_arr


res0Arr = residualscalc(file_path0)
res1Arr = residualscalc(file_path1)
print("LAPPDS 0 RESIDUALS", res0Arr)
print("LAPPDS 1 RESIDUALS ", res1Arr)
ak.flatten(res0Arr)
ak.flatten(res1Arr)
both_flat = ak.concatenate([ak.flatten(res0Arr), ak.flatten(res1Arr)])
print(len(both_flat))


plt.hist(both_flat, bins = 100, range = (-80, 20))
plt.xlabel("time(ns)")
plt.ylabel("count")
plt.title("All Residuals n over thresh > 1 and amp > 10V")
plt.show()






    


