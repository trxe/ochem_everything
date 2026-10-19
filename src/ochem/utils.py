import numpy as np
import os
import pandas as pd
from matplotlib import pyplot as plt
from pathlib import Path
import statistics as stat

def get_samples(data,
                expno, #used to generate list of sample names. do not use "_"!
                rxns, #list of rxns you run in this expt
                n_rows = 6,
                n_cols = 8,
                export_list = True,
                export_map = False):
      
    if len(data) > 2 * n_cols:
        df = data.loc[0:2 * n_cols - 1]
        print("Too many time points! Sample list cut off from "+str(data["time"][2*n_cols])+"m onwards!")
    else:
        df = data
        
    timepts = df["time"]
      
    platerows = [chr(65+i) for i in range(n_rows)]

    vialmap = pd.DataFrame()
    namemap = pd.DataFrame()
    vials = []
    samplenames = []

    for rxn in rxns:
        vialmap[rxn] = []
        namemap[rxn] = []

    for i in range(len(timepts)):
        
        vialrow = []
        namerow = []
        
        for j in range(len(rxns)):
            vialname = ("1:" if i < n_cols else "2:") + platerows[j] + "," + str(i % n_cols + 1)
            vials.append(vialname)
            vialrow.append(vialname)
            
            samplename = "_".join([expno,rxns[j],str(timepts[i])]) #format: EXPNO_RXN_TIME
            samplenames.append(samplename)
            namerow.append(samplename)
            
        vialmap.loc[len(vialmap)] = vialrow
        namemap.loc[len(namemap)] = namerow
    
    if export_list:        
        with open("-".join([expno,"vialnames.txt"]),"w") as file:
            for vial in vials:
                file.write(vial+"\n")
                
        with open("-".join([expno,"samplenames.txt"]),"w") as file:
            for samplename in samplenames:
                file.write(samplename+"\n")
    
    if export_map:
        vialmap["time"] = timepts
        vialmap.to_csv(expno+"-vial-map.csv", index = False)
    
        namemap['time'] = timepts
        namemap['temp'] = data['temp']
        namemap.to_csv(expno+"-name-map.csv", index = False)
    
'''
user df structure:
name - RT - tolerance

waters df structure
"RT" "SampleName" "Area" "Base Peak"

'''
def extract_peaks(export_data, substance_list, name_map, export = "default", channel = "UV"): #todo: extract sample names from name map
    
    rxnlist = name_map.keys().tolist() # get list of rxns from name map
    rxnlist.remove("time")
    rxnlist.remove("temp")
    #rxnlist.remove('Unnamed: 0') # get rid of non-rxn name columns
    
    timepts = name_map["time"].tolist()
    temppts = name_map["temp"].tolist() #save time and temp columns of name map for later use
    
    for i in range(len(substance_list)):
        
        substance_peak_list = pd.DataFrame({ #initialize df for each substance with time and temp
            "time": timepts,
            "temp": temppts
        })
        
        RT_min = substance_list["RT"][i] - substance_list["tolerance"][i],
        RT_max = substance_list["RT"][i] + substance_list["tolerance"][i]
        
        candidate_peaks = pd.DataFrame()
        
        samplenames = []
        
        for rxn in rxnlist:
            substance_peak_list[rxn] = np.nan #initialize columns of df with rxn names
            samplenames += name_map[rxn].tolist() #get a long 1 dimensional list of all sample names from the name map
    
        candidate_peaks = export_data[(export_data["RT"] >= RT_min) & (export_data["RT"] <= RT_max) & (export_data["Channel Name"] == channel)] #goes through export data, picks out every candidate peak with RT matching criteria
        
        candidate_peaks = candidate_peaks.reset_index(drop = True) #make candidate_peaks indices start from 0 again
        
        for samplename in samplenames: #iterate thru all samples to find peaks for each one
            
            expno, rxn, time = samplename.split("_")
            row_index = timepts.index(int(time))
            
            possible_peaks = []
            
            for k in range(len(candidate_peaks)):
                
                if candidate_peaks["SampleName"][k] == samplename: #now list all the peak areas with SampleName that matches
                    
                    possible_peaks.append(candidate_peaks["Area"][k])
            
            if len(possible_peaks) > 0:
                substance_peak_list[rxn][row_index] = max(possible_peaks) #select biggest peak and write it into the peak list
            else:00000
pvis5761.40.374308489785250660.2263157894736842200.144333725320294620
                substance_peak_list[rxn][row_index] = 0 #if no possible peaks detected for this sample, set peak area to 0
                    
        if export == "default":
            substance_peak_list.to_csv("_".join([expno,substance_list["name"][i]])+".csv") #write a new csv with name format EXPNO_SUBSTANCE.csv for each substance detected
        elif export:
            substance_peak_list.to_csv("_".join([export,substance_list["name"][i]])+".csv")
            
def plot_ramp(data, rxns, #rxn names from columns in spreadsheet
              y_name = None, #y-axis label
              y_labels = None, #plot labels, same as reaction names by default
              y_range = None, #y-axis scale 
              ramp = True, #True if not isothermal (has temp data)
              temp_range = None, #temperature range, x-axis scale
              colors = ["red","darkorange","gold","mediumseagreen","royalblue","darkviolet"], #default line colors
              title = ""
              ):
    
    fig, ax1 = plt.subplots()
    
    time = data["time"]    
    
    ax1.set_xlabel("Time", weight="bold")
    
    if not y_labels:
        y_labels = rxns
    
    '''
    if not y_range: #find max and min y values in selected columns
        maxls = []
        minls = []
        
        for label in y_labels:
            maxls.append(max(data[label]))
            minls.append(min(data[label]))
            
        max_y = max(maxls)
        min_y = min(minls)
        
        y_range = [min_y - (max_y - min_y) * 0.1,
                   max_y + (max_y - min_y) * 0.1] #sets default y_range to range of selected columns +10% on either end
    '''
    
    if y_range:
        ax1.set_ylim(y_range[0],y_range[1]) #y-axis data scale
    
    for i in range(len(rxns)):
        ax1.plot(time, data[rxns[i]], colors[i], label=y_labels[i], marker = "o")
        
    if y_name:
        ax1.set_ylabel(y_name, weight="bold")
    
    plt.title(title)
    plt.legend(loc='upper left', bbox_to_anchor=(1.15, 1.0)) #change location of legend here
    
    if ramp:
        temp = data["temp"]
        ax2 = ax1.twinx()
        ax2.set_ylabel("Temperature", weight="bold")
        
        if not temp_range:
            temp_range = [min(temp)-5,
                          max(temp)+5]
        
        ax2.set_ylim(temp_range[0],temp_range[1]) #y-axis temperature scale
        
        ax2.plot(time,temp,"grey")

    plt.show()

def rate_from_ramp(data, y_labels, export: Optional[str] = None):
    data['delta_time'] = data['time'] - data['time'].shift(1)
    if 'temp' in data.columns:
        data['avg_temp'] = (data['temp'] - data['temp'].shift(1)) / 2
    for label in y_labels:
        data[f"rate_{label}"] = (data[label] - data[label].shift(1)) / data['delta_time']
    if export:
        data.to_csv()
    return data

def area_to_conc(data, y_labels, by, export = False):
    
    df = pd.DataFrame()
    
    df["temp"] = data["temp"]
    
    df["time"] = data["time"]
    
    for label in y_labels:
        df[label] = data[label] * by
        
    if export:
        df.to_csv(export)
        
    return df

def a_by_b(data, y_labels, by, export = False):
    
    df = pd.DataFrame()
    
    df["temp"] = data["temp"]
    
    df["time"] = data["time"]
    
    for label in y_labels:
        
        df[label] = data[label] / by[label]
        
    if export:
        
        df.to_csv(export)
        
    return df

def moving_avg(data, y_labels, window, export = False):
    cols = ['time', 'temp', *y_labels]
    df = data.copy()

    for col in cols:
        df[col] = sum([df[col].shift(i) for i in range(window)]) / window

    if export:
        df.to_csv(export)
    
    return df

def VTNA(data, #data containing (usually) pdt concentrations
         rxns,
         concs, #list of dfs of concentrations of reactant trying to do VTNA on
         powers, #assumes first order in reaction, pass as list e.g. [1, 0.5, 2, 3], must same length as concs,
         y_labels = None,
         title = None,
         y_name = None,
         x_name = None, #if you want a x label, provide a list of reagent names in the same order as concs and powers above here
         plot = True,
         ramp = True): #doesn't do anything, temperature axis is meaningless in VTNA
    
    time = data["time"]
    
    df = pd.DataFrame({
        "time":data["time"],
    })
    
    if ramp:
        df["temp"] = data["temp"]
    
    ls = []

    for rxn in rxns:
        ls.append([0])
    
    for i in range(1,len(time)): #iterates thru time points
        
        deltaT = time[i] - time[i-1]
        
        for j in range(len(rxns)): #at each time pt, iterate thru rxns
            
            summation = ls[j][i-1] # starting value of value to be added
            product = 1
            
            for k in range(len(concs)):
                mean_conc = (concs[k][rxns[j]][i] + concs[k][rxns[j]][i-1]) / 2
                
                product *= np.power(mean_conc, float(powers[k]))

            ls[j].append(summation+product*deltaT)
    
    if plot:
        fig, ax = plt.subplots(figsize = (8,4))
        
        for i in range(len(rxns)):
            plt.plot(ls[i],data[rxns[i]],label = rxns[i] if not y_labels else y_labels[i], marker = "o")
        
        if title:
            ax.set_title(title, weight = "bold")
        
        if y_name:
            ax.set_ylabel(y_name, weight = "bold")
            
        if x_name:
            x_label = "\u2211"
            table = str.maketrans("-.0123456789", "⁻˙⁰¹²³⁴⁵⁶⁷⁸⁹")
            for i in range(len(x_name)):
                x_label = x_label + "[" + x_name[i] + "]" + str(powers[i]).translate(table)
            x_label += "\u0394t"
            ax.set_xlabel(x_label, weight = "bold", fontsize = 12)
            
        plt.legend()
        
        plt.show()
    
    for i in range(len(rxns)): #convert list to df for export
        df[rxns[i]] = ls[i]
    
    return df

'''
ALL UNITS IN M BEFORE PASSING TO THIS FUNCTION!
'''
def eyring(data, rxns,
           T = 298, #default temperature to give deltaG
           points = False, #uses all data points by default
           export = False
           ):
    
    df = pd.DataFrame()
    
    kB = 1.381 * (10 ** -23)
    h = 6.626 * (10 ** -34)
    R = 1.987
    
    if points:
        data = data.iloc[points[0]:points[1],]
    
    df["temp"] = data["temp"] + 273.15 #converts all temperature to K assuming C
    df["invtemp"] = 1 / (df["temp"])
    templs = df["invtemp"].tolist()

    fig, ax = plt.subplots()

    for rxn in rxns:
        
        df[rxn] = np.log((data[rxn] / 60) / df["temp"]) #assumes rates in data are in min^-1, takes ln(k/T) in seconds
                             
        plt.scatter(templs, df[rxn], label=rxn)

    plt.ylabel("ln(k/T)",weight="bold")
    plt.xlabel("1/T",weight="bold")
    ax.ticklabel_format(useMathText=True,scilimits=[-2,2])

    datalist = [[],[]]
    
    for rxn in rxns:
        datalist[0] += templs
        datalist[1] += df[rxn].tolist()

    idx = np.isfinite(datalist[0]) & np.isfinite(datalist[1])
    
    long_df = pd.DataFrame({"invtemp":datalist[0],
                            "y":datalist[1]})
    
    temp_all = long_df["invtemp"][idx]
    y_all = long_df["y"][idx]
    
    linearfit = np.polyfit(temp_all, y_all,1)
    print(linearfit)
    line = np.poly1d(linearfit)
    print(line)

    deltaH = -linearfit[0] * R / 1000
    print("\u0394H\u2021 = " + str(deltaH))

    deltaS = (linearfit[1] - np.log(kB/h)) * R
    print("\u0394S\u2021 = "+ str(deltaS))

    deltaG = deltaH - T * deltaS / 1000
    print("\u0394G\u2021 = " + str(deltaG))

    x_line = np.linspace(min(temp_all), max(temp_all), 100)
    plt.plot(x_line,line(x_line))

    y_pred = line(temp_all)
    residuals = y_all - y_pred
    ss_res = sum(residuals ** 2)
    ss_tot = sum(((y_all-np.mean(y_all))**2))
    r_sq = 1 - (ss_res / ss_tot)

    print(r_sq)

    plt.annotate("y = " + str(round(linearfit[0],3)) + "x + " + str(round(linearfit[1],3)) +
                "\nR^2 = "+str(round(r_sq,3))+
                "\n\u0394H\u2021 = " + str(round(deltaH,3)) + " kcal/mol\n"+
                "\u0394S\u2021 = " + str(round(deltaS,3)) + " cal/mol\n"+
                "\u0394G\u2021" + "(" + str(T) + "K) = " + str(round(deltaG,3)) + " kcal/mol",
                xy=(0.05, 0.075),
                xycoords="axes fraction",
                fontsize=10,
                bbox=dict(boxstyle="round,pad=0.5",facecolor = "white", edgecolor="gray")
    )

    #ax.set_xticklabels([])
    plt.title("Eyring plot for product", weight="bold")

    plt.legend(loc = "upper right", bbox_to_anchor = (0.975,0.975), edgecolor="gray")
    plt.show()

    plt.scatter(temp_all,residuals)
    
    plt.plot(x_line,np.poly1d([0,0])(x_line),color="black",linestyle=(0, (1, 10)))
    plt.title("Residuals",weight="bold")
    plt.ticklabel_format(useMathText=True,scilimits=[-2,2])
    plt.xlabel("1/T",weight="bold")
    plt.show()
    
    if export:
        df.to_csv(export)

def convert_autoVTNA(datas, #list of dataframes with output data
                     rxnnames, #rxn names as listed in data
                     substnames, #list of substrate names to be passed to auto-VTNA, in the same order as datas above
                     export = False
                     ):
    
    maindict = {}
    
    for rxnname in rxnnames:
        
        rxndict = {}
        
        for i in range(len(datas)):
            
            rxndict["time"] = datas[i]["time"]
            rxndict[substnames[i]] = datas[i][rxnname]
        
        if export:
            
            df = pd.DataFrame(rxndict)
            df.to_csv(export+"_"+rxnname+".csv", index=False)

    
