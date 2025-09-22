import uproot
import awkward as ak 
import numpy as np
import matplotlib.pyplot as plt

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
    center = []
    raw = []
    trigs = []

    
    if "lappd_0" in file_path:
        for i in range(len(waveforms)):
            evtgroup0 = []
            evtgroup1 = []
            evtgroup2 = []
            evtgroup3 = []
            evttrig = []
            center_event07 =[]
            center_event8 = []


            for j in range(len(waveforms[i])):
                if lcns[i][j] < 31:
                    ped = ak.mean(waveforms[i][j][500:900])
                    wf = waveforms[i][j][10:-10]
                    amp = ak.max(wf) - ped
                    if amp > 75:
                        crossed = wf - ped > (0.6 * (ak.max(wf) - ped))
                        above = ak.where(crossed)[0]
                        if len(above) > 0:
                            xpt_above = above[0]
                            xpt_below = xpt_above - 1
                            y_below = wf[xpt_below]
                            y_above = wf[xpt_above]
                            threshold = 0.6 * (ak.max(waveforms[i][j][10:-10]) - ped) + ped
                            linearized_index = linear_interpolation(xpt_below, y_below, xpt_above, y_above, threshold)
                            linearized_index_ns =  (linearized_index + 10) * 0.2 ## account for shift 
                            if lcns[i][j] not in [0, 7, 8]:  
                                if lcns[i][j] < 8:
                                    evtgroup0.append(linearized_index_ns)
                                if 7 < lcns[i][j] < 16:
                                    evtgroup1.append(linearized_index_ns)
                                if 15 < lcns[i][j] < 24:
                                    evtgroup2.append(linearized_index_ns)
                                if 23 < lcns[i][j] < 31:  # 31 is coincident PMT trigger
                                    evtgroup3.append(linearized_index_ns)
                            if lcns[i][j]  in [0, 7]:
                                center_event07.append(linearized_index_ns)
                            if lcns[i][j] == 8: 
                                center_event8.append(linearized_index_ns)



                elif lcns[i][j] > 31:  # trigger PMTs 32, 48, 64, 80
                    ped = ak.mean(waveforms[i][j][100:500])
                    wf = waveforms[i][j][100:-100]
                    amp = ak.max(wf - ped)
                    if amp > 50:
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
                    raw.append(t)
                for t in evtgroup1:
                    tres.append(t - evttrig[1])
                    raw.append(t)
                for t in evtgroup2:
                    tres.append(t - evttrig[2])
                    raw.append(t)
                for t in evtgroup3:
                    tres.append(t - evttrig[3])
                    raw.append(t)
                for t in center_event07: 
                    center.append(t - evttrig[0])
                for t in center_event8: 
                    center.append(t - evttrig[1])

        return [tres, center]
    
    
    
    
    
    if "lappd_1" in file_path:
        for i in range(len(waveforms)):
            evtgroup0 = []
            evtgroup1 = []
            evtgroup2 = []
            evtgroup3 = []
            evttrig = []
            center_event = []

            for j in range(len(waveforms[i])):
                if lcns[i][j] < 31:
                    ped = ak.mean(waveforms[i][j][500:900])
                    wf = waveforms[i][j][10:-10]
                    amp = ak.max(wf) - ped
                    if amp > 75:
                        crossed = wf - ped > (0.6 * (ak.max(wf) - ped))
                        above = ak.where(crossed)[0]
                        if len(above) > 0:
                            xpt_above = above[0]
                            xpt_below = xpt_above - 1
                            y_below = wf[xpt_below]
                            y_above = wf[xpt_above]
                            threshold = 0.6 * (ak.max(waveforms[i][j][10:-10]) - ped) + ped
                            linearized_index = linear_interpolation(xpt_below, y_below, xpt_above, y_above, threshold)
                            linearized_index_ns =  (linearized_index + 10) * 0.2 ## account for shift 
                            if lcns[i][j] != 25:
                                if lcns[i][j] < 8:
                                    evtgroup0.append(linearized_index_ns)
                                if 7 < lcns[i][j] < 16:
                                    evtgroup1.append(linearized_index_ns)
                                if 15 < lcns[i][j] < 24:
                                    evtgroup2.append(linearized_index_ns)
                                if 23 < lcns[i][j] < 31:  # 31 is coincident PMT trigger
                                    evtgroup3.append(linearized_index_ns)

                            if lcns[i][j] == 25:
                                center_event.append(linearized_index_ns)
                
                
                elif lcns[i][j] > 31:  # trigger PMTs 32, 48, 64, 80
                    ped = ak.mean(waveforms[i][j][100:500])
                    wf = waveforms[i][j][100:-100]
                    amp = ak.max(wf - ped)
                    if amp > 50:
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
                    raw.append(t)
                for t in evtgroup1:
                    tres.append(t - evttrig[1])
                    raw.append(t)
                for t in evtgroup2:
                    tres.append(t - evttrig[2])
                    raw.append(t)
                for t in evtgroup3:
                    tres.append(t - evttrig[3])
                    raw.append(t)
                for t in center_event:
                    center.append(t - evttrig[3])
        return [tres, center]










residuals_scint_scintillator = residualcalcman(paths_scint[0])[0] + residualcalcman(paths_scint[1])[0] + residualcalcman(paths_scint[2])[0] + residualcalcman(paths_scint[3])[0] 
residuals_scint_cherenkov = residualcalcman(paths_scint[0])[1] + residualcalcman(paths_scint[1])[1] + residualcalcman(paths_scint[2])[1] + residualcalcman(paths_scint[3])[1]
# Plot residual histograms for water and scint
plt.figure(figsize=(8,6))

bins = 200  
range_min, range_max = -100, 10  #

plt.hist(residuals_scint_scintillator, bins=bins, range=(range_min, range_max),
         histtype='step', color='blue', label='scintillator', linewidth=1.5)

plt.hist(residuals_scint_cherenkov, bins=bins, range=(range_min, range_max),
         histtype='step', color='red', label='cherenkov', linewidth=1.5)

plt.xlabel("Residual [ns]")
plt.ylabel("Counts")
plt.legend()
plt.show()











