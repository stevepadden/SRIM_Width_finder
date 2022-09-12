#imports
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from gaussians import *
from scipy import interpolate

#Setting up imports
folder = r"D:\Steve\Documents\SRIM_OUTPUT\\"    #Srim output files
prefix = r"SRIM_SI_100_5_"  #File prefix WITHOUT thickness number
suffix = r".txt"    #File suffix
num_sim = 5000 # Number of simulated particles - keep at 5k in both this and srim for ease

ranges = [6945,7445,7945,8445,8945,9445,9945,10445,10945,8500,8600,8700,8800,8571]  #These are both the ranges and the Thickness number appended to the file path!

ranges.sort()
#column names
transcols = ["Ion Numb","Atom Numb","Energy","Depth_X","Lateral_Y","Lateral_Z","Cos_X","Cos_Y","Cos_Z"]
#Function that finds the percentages for each SRIM file
def Find_Percentages(range):
    data = (folder+prefix+"%s"+suffix) % range
    with open(data, 'r') as file:
        text = file.read()
    clean = text.replace("T", "")

    with open(folder + "cleaned.txt", "w+") as file:
        file.write(clean)
    try:
        data = np.genfromtxt(folder + "cleaned.txt", skip_header=12, invalid_raise=False, dtype="float")
        print("INCOMING DATA FOR RANGE %s" % range)
        data = pd.DataFrame(data, columns=transcols)
        print(data)
        tp = (len(data["Energy"]) / num_sim) * 100  # The percentage of particles transmitted through the foil

        sube = data["Energy"][data["Energy"] < 5000]
        numsube = sube.count()
        sub5p = (numsube / num_sim) * 100


        #axhist.hist(data["Energy"],bins=50,edgecolor=(0,0,0,0.4),label="Silicon:%s" % (range))

        return tp, sub5p
    except:
        print("No Data found in range file %s ; Returning 0" % (range))
        return 0,0

def plot_hists(range,alpha=1):
    data = (folder+prefix+"%s"+suffix) % range
    with open(data, 'r') as file:
        text = file.read()
    clean = text.replace("T", "")

    with open(folder + "cleaned.txt", "w+") as file:
        file.write(clean)
    try:
        data = np.genfromtxt(folder + "cleaned.txt", skip_header=12, invalid_raise=False, dtype="float")
        print("INCOMING DATA FOR RANGE %s" % range)
        data = pd.DataFrame(data, columns=transcols)
        axhist.hist(data["Energy"].div(1000),bins=50,edgecolor=(0,0,0,0.4),alpha=alpha,label=("Silicon:%s $\AA$" % (range)))
    except:
        print("No Data found in range file %s ; Returning 0" % (range))
        return 0, 0

def overlap_hists(range,alpha=1,switch=1):
    data = (folder+prefix+"%s"+suffix) % range
    with open(data, 'r') as file:
        text = file.read()
    clean = text.replace("T", "")

    with open(folder + "cleaned.txt", "w+") as file:
        file.write(clean)
    try:
        data = np.genfromtxt(folder + "cleaned.txt", skip_header=12, invalid_raise=False, dtype="float")
        print("INCOMING DATA FOR RANGE %s" % range)
        data = pd.DataFrame(data, columns=transcols)
        if switch == 1:
            axhist.hist(data["Energy"].div(1000),bins=50,histtype='step',edgecolor=(0,0,0,1),alpha=alpha,label=("Silicon:%s $\AA$" % (range)),lw=1.3)
        else:
            axhist.hist(data["Energy"].div(1000), bins=50, histtype='step', edgecolor=(1,1,1, 1), alpha=alpha,lw=1.3,
                        label=("Silicon:%s $\AA$" % (range))) # (0.7,0.7,0.7, 1)
    except:
        print("No Data found in range file %s ; Returning 0" % (range))
        return 0, 0

Total_Percent = []
Sub5_Percent = []

for i in ranges:
    total, sub= Find_Percentages(i)
    Total_Percent.append(total)
    Sub5_Percent.append(sub)




submodel = Model(gaussian)  # define the function returned model
Amp = 50
Cen = 8571
Wid = 683
xnew = np.arange(min(ranges),max(ranges),50)
result = submodel.fit(Sub5_Percent, x=ranges, amp=Amp, cen=Cen, wid=Wid,method="least_squares")
tamp = result.params['amp'].value
twid = result.params['wid'].value
tcen = result.params['cen'].value
const_fit = gaussian(xnew, tamp, tcen, twid)  # y values from model
print("Suggested Width : %f"%tcen)
print(result.fit_report())

totmodel = Model(inverse_s)
head =100
a = 9000
b=5
tck = interpolate.splrep(ranges,Total_Percent)
spline = interpolate.splev(xnew,tck)





fit = mgauss(xnew,Amp,Cen,Wid)

fig,ax = plt.subplots(figsize=(6.25,4))

t1 =ax.scatter(ranges,Total_Percent,marker="+",color="royalblue")
s1 =ax.scatter(ranges,Sub5_Percent,marker="+",color="red")
s2=ax.plot(xnew,const_fit,color="red",label="Percentage transmitted \nwith energy $< 5$keV")
#ax.plot(xnew,const_tot_fit,color="blue")
t2=ax.plot(xnew,spline,label="Total percentage \ntransmitted through foil",linestyle="--",color="royalblue")
#ax.plot(ranges,result.best_fit)
#ax.plot(xnew,fit)
#ax.plot(ranges,result.init_fit)
ax.set_xlabel("Foil Thickness ($\AA$)")
ax.set_ylabel (r"Percentage (%)")
ax.set_xlim((min(ranges) ,max(ranges)))
ax.set_title("SRIM 100 keV protons into silicon")
str = "Optimal thickness: %i $\AA$" % (8571)#tcen

print(ranges)
ax.axvline(8570,color="#878b91",alpha=0.4,label=str)
ax.scatter(ranges[5],Sub5_Percent[5],marker="^",color="red",label="Optimal thickness simulation")
ax.legend(loc="upper right",frameon=False)#,fontsize=10)
#ax.legend([(t1,t2),(s1,s2)],['set1','set2'])
# import matplotlib.lines as mlines
# orange_line = mlines.Line2D([], [], color='orange', marker='+',
#                           markersize=10, label='Percentage transmitted with energy $< 5$keV')
# ax.legend(handles=[orange_line])
fig.show()
df = pd.DataFrame({"ranges":ranges,"Total":Total_Percent,"Sub5":Sub5_Percent})
df.to_csv("Transmissions.csv",index=False)


reduced_ranges = [6945,7445,7945,8445,8571]
reduced_ranges.sort()#(reverse=True)
fighist,axhist = plt.subplots()
axhist.set_xlabel("Energy of transmitted protons (keV)")
axhist.set_ylabel("Counts")
axhist.axvline(5,color="#878b91",alpha=1,label="5keV")
axhist.set_xlim((0,30))
axhist.set_title("SRIM simulated energy distributions of transmitted protons")
for i in reduced_ranges:
    plot_hists(i,alpha=0.6)
axhist.legend(loc="upper right",frameon=False)

reduced_ranges.sort(reverse=True)
count = 0
for i in reduced_ranges:
    if count%2 == 0:
        switch = 0
        print(switch)
    else:
        switch = 1
        print(switch)
    overlap_hists(i,1,switch)
    count = count+1



fighist.show()
print(ranges)
print(Sub5_Percent)