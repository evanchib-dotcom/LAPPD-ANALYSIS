import uproot
import awkward as ak 
import numpy as np
import matplotlib.pyplot as plt

file_path1 = "C:\\Users\\Evan\\Downloads\\lappd-158-water-beta-mask-500nm-lp-centered_lappd_1.root"
file_path0 = "C:\\Users\\Evan\\Downloads\\lappd-158-water-beta-mask-500nm-lp-centered_lappd_0.root"




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


    n = 1000
    tres = []
    raw = []
    trigs = []

    for i in range(len(waveforms[:n])):
        evtgroup0 = []
        evtgroup1 = []
        evtgroup2 = []
        evtgroup3 = []
        evttrig = []

        for j in range(len(waveforms[i])):
            if lcns[i][j] < 31:
                ped = ak.mean(waveforms[i][j][500:900])
                wf = waveforms[i][j][10:-10]
                amp = ak.max(wf - ped)
                if amp > 50:
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
    return tres



amp_mask = 10   ## in Volts
## residual calculating usingg already processed times
def residualscalc(file_path):
    n_tup = uproot.open(file_path)
    output = n_tup['output;1']
    event_info = output.arrays(output.keys(),library = "ak")
    waveforms = event_info['waveforms']
    peaks = event_info['peaks'] 
    ped = event_info['pedestals'] 
    lcns =  event_info['lcns']
    n_abv_thresh = event_info['nsamples_above_threshold']
    times = event_info['times']
    t_res = []
    discard = []
    center_pixels = []
    

    
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
        
        
    if file_path == file_path0 :   
        for i in range(len(waveforms)):
            discard_channel_per_event = []
            center_pixels_evt = []
            t_res_event = []
            for j in range(len(waveforms[i])):
                amp = peaks[i][j] - ped[i][j]
                if (n_abv_thresh[i][j] > 1) and (amp > amp_mask) :
                    if lcns[i][j] not in [31,32,48,64,80]: ## not the trigger or PMT 
                        residual = times[i][j] - locate_groputrigger_time(i,j) 
                        if lcns[i][j] not in [0,7,8]: ## not the center pixels 
                            t_res_event.append(residual)
                        elif lcns[i][j] in [0,7,8]: ## if they are center pixels append them to different list
                            center_pixels_evt.append(residual)
                    
                    else:
                        discard_channel_per_event.append(j)
                else : 
                    discard_channel_per_event.append(j)

            discard.append(discard_channel_per_event)
            t_res.append(t_res_event) 
            center_pixels.append(center_pixels_evt)

      
       
    
    
    
    if file_path == file_path1 :
         for i in range(len(waveforms)):
            discard_channel_per_event = []
            center_pixels_evt = []
            t_res_event = []
            for j in range(len(waveforms[i])):
                amp = peaks[i][j] - ped[i][j]
                if (n_abv_thresh[i][j] > 1) and (amp > amp_mask) :
                    if lcns[i][j] not in [31,32,48,64,80]: ## not the trigger or PMT 
                        residual = times[i][j] - locate_groputrigger_time(i,j) 
                        if lcns[i][j] != 25: ## not the center pixels 
                            t_res_event.append(residual)
                        elif lcns[i][j] == 25: ## if they are center pixels append them to different list
                            center_pixels_evt.append(residual)
                    
                    else:
                        discard_channel_per_event.append(j)
                else : 
                    discard_channel_per_event.append(j)

            discard.append(discard_channel_per_event)
            t_res.append(t_res_event) 
            center_pixels.append(center_pixels_evt)

    t_res_arr = ak.Array(t_res)
    center_arr = ak.Array(center_pixels)    
    return [t_res_arr, center_arr]
    







res_without_center1 = residualscalc(file_path1)[0]
res_without_center0 = residualscalc(file_path0)[0]

both_lappds_nocenter = ak.concatenate([ak.flatten(res_without_center0),ak.flatten(res_without_center1)])

res_center0 = residualscalc(file_path0)[1]
res_center1 = residualscalc(file_path1)[1]

both_lappds_center = ak.concatenate([ak.flatten(res_center0),ak.flatten(res_center1)])

print("long wavelength pixels residuals", both_lappds_nocenter)
print("# of residuals" , len(both_lappds_nocenter))

print("short wavelength pixels residuals", both_lappds_center)
print("# of residuals" , len(both_lappds_center))



plt.hist(both_lappds_center, bins = 100, range = (-80, 20))
plt.xlabel("time(ns)")
plt.ylabel("count")
plt.title("Residuals center pixels only 158 mask, n > 1 and amp > 10V")
plt.show()



plt.hist(both_lappds_nocenter, bins = 100, range = (-80, 20))
plt.xlabel("time(ns)")
plt.ylabel("count")
plt.title("Residuals everything but center pixels 158,  mask n > 1 and amp > 10V")
plt.show()
