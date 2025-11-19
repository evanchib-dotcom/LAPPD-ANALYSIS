import uproot
import awkward as ak 
import numpy as np
import matplotlib.pyplot as plt



paths_lab = [
    "C:\\Users\\Evan\\Downloads\\lappd-158-lab1.7ppo-beta-no-mask-band-blocker-500nm-lp-centered-redo_lappd_1.root",
    "C:\\Users\\Evan\\Downloads\\lappd-158-lab1.7ppo-beta-mask-band-blocker-500nm-lp-centered_lappd_0.root",
    "C:\\Users\\Evan\\Downloads\\lappd-158-lab1.7ppo-beta-mask-band-blocker-500nm-lp-centered_lappd_1.root",
    "C:\\Users\\Evan\\Downloads\\lappd-158-lab1.7ppo-beta-no-mask-band-blocker-500nm-lp-centered-redo_lappd_0.root",]
    
paths_water = [    
     "C:\\Users\\Evan\\Downloads\\lappd-158-water-beta-no-mask-band-blocker-500nm-lp-centered_lappd_1.root",
    "C:\\Users\\Evan\\Downloads\\lappd-158-water-beta-no-mask-band-blocker-500nm-lp-centered_lappd_0.root",
    "C:\\Users\\Evan\\Downloads\\lappd-158-water-beta-no-mask-500nm-lp-centered_lappd_0.root",
    "C:\\Users\\Evan\\Downloads\\lappd-158-water-beta-no-mask-500nm-lp-centered_lappd_1.root",
    "C:\\Users\\Evan\\Downloads\\lappd-158-water-beta-mask-500nm-lp-centered_lappd_0.root",
    "C:\\Users\\Evan\\Downloads\\lappd-158-water-beta-mask-500nm-lp-centered_lappd_1.root", ]

path_alt = "C:/Users/Evan/Downloads/lappd-lab-beta-no-mask-band-blocker-500nm-lp-centered-redo-125-lappd_1.root"
path_alt0 = "C:\\Users\\Evan\\Downloads\\lappd-lab-beta-no-mask-band-blocker-500nm-lp-centered-redo-125-lappd_0.root"
amps = []

def residualcalcman(file_path):
    n_tup = uproot.open(file_path)
    output = n_tup["output;1"]
    event_info = output.arrays(output.keys(), library="ak")
    lcns = event_info['lcns']
    waveforms = event_info['waveforms']


    def linear_interpolation(xbel, ybel, xab, yab, threshold):
        
        if yab == ybel:
            return xbel  # Avoid division by zero, just return first point
        return xbel + (threshold - ybel) * (xab - xbel) / (yab - ybel)


   
    tres = []
    trigs = []
    for i in range(len(waveforms)):
            evtgroup0 = []
            evtgroup1 = []
            evtgroup2 = []
            evtgroup3 = []
            evttrig = []
            


            for j in range(len(waveforms[i])):
                if lcns[i][j] < 31:
                    ped = ak.mean(waveforms[i][j][500:900])
                    wf = waveforms[i][j][10:-10]
                    amp = ak.max(wf) - ped
                    print(amp)
                    if amp > 55:
                        crossed = wf - ped > (0.6 * (ak.max(wf) - ped))
                        above = ak.where(crossed)[0]
                        if len(above) > 0:
                            xpt_above = above[0]
                            xpt_below = xpt_above - 1
                            y_below = wf[xpt_below]
                            y_above = wf[xpt_above]
                            threshold = 0.6 * (ak.max(waveforms[i][j][10:-10]) - ped) + ped
                            linearized_index = linear_interpolation(xpt_below, y_below, xpt_above, y_above, threshold)
                            linearized_index_ns =  (linearized_index + 10) * 0.2 
                    
                            if lcns[i][j] < 8:
                                    evtgroup0.append(linearized_index_ns)
                            if 7 < lcns[i][j] < 16:
                                    evtgroup1.append(linearized_index_ns)
                            if 15 < lcns[i][j] < 24:
                                    evtgroup2.append(linearized_index_ns)
                            if 23 < lcns[i][j] < 31:  # 31 is coincident PMT trigger
                                    evtgroup3.append(linearized_index_ns)
                          



                elif lcns[i][j] > 31:  # trigger PMTs 32, 48, 64, 80
                    ped = ak.mean(waveforms[i][j][100:500])
                    wf = waveforms[i][j][100:-100]
                    amp = ak.max(wf - ped)
                    crossed = wf - ped > (0.6 * (ak.max(wf) - ped))
                    above = ak.where(crossed)[0]
                    if len(above) > 0:
                            threshold = 0.6 * (ak.max(wf) - ped) + ped 
                            linearized_trigger_t = linear_interpolation(above[0] - 1, wf[above[0] - 1], above[0],wf[above[0]] , threshold)
                            evttrig.append(0.2 * (linearized_trigger_t + 100))
                          
                
           
           
            if len(evttrig) == 4:
                for t in evttrig:
                    trigs.append(t)
                for t in evtgroup0:
                    tres.append(t - evttrig[0])
                for t in evtgroup1:
                    tres.append(t - evttrig[1])
                for t in evtgroup2:
                    tres.append(t - evttrig[2])
                for t in evtgroup3:
                    tres.append(t - evttrig[3])   

          
    return tres







n_tup = uproot.open(path_alt0)
output = n_tup["output;1"]
event_info = output.arrays(output.keys(), library="ak")
lcns = event_info['lcns']
waveforms = event_info['waveforms']



residuals_lab = residualcalcman(path_alt)


plt.figure(figsize=(8,6))





bins = 250
range_min, range_max = -250, 250



plt.hist(residuals_lab, bins=bins, range=(range_min, range_max),
histtype='step')

plt.xlabel(" time [ns]")
plt.ylabel("Counts")
plt.title("python residuls dig1")
plt.show()
