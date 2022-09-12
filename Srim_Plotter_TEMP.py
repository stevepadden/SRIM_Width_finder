#Imports
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt



#Importing the transmitted ions and huge depth
folder = r"D:\Steve\Documents\SRIM_OUTPUT\\"
depth = np.genfromtxt(folder+"RANGE_3D_Si.txt",skip_header=17)
transmitted_string = folder+"TRANSMIT_100_5_8571_SI.txt"
#transmitted_string = folder+"TRANSMIT_100_5_7945_SI.txt"
total_num_sim = 30000
Foil_Width=8571#
PRAL_thickness=7945
#function needed for square images
def findmax(x, y):
    x1 = max(x)
    x2 = abs(min(x))
    y1 = max(y)
    y2 = abs(min(y))
    if x2 > x1:
        maxx = x2
        minx = -x2
    else:
        maxx = x1
        minx = -x1

    if y2 > y1:
        maxy = y2
        miny = -y2
    else:
        maxy = y1
        miny = -y1
    endlim = (max(maxx, maxy))
    return minx, maxx, miny, maxy, endlim

#Trimming the weirdly broken TRIM output files
with open(transmitted_string,'r') as file:
    text = file.read()
clean = text.replace("T","")

with open(folder+"cleaned.txt","w+") as file:
    file.write(clean)

transmitted = np.genfromtxt(folder+"cleaned.txt",skip_header=12,invalid_raise=False,dtype="float")
depthcols = ["Ion_Number","Depth_X","Lateral_Y","Lateral_Z"]
transcols = ["Ion Numb","Atom Numb","Energy","Depth_X","Lateral_Y","Lateral_Z","Cos_X","Cos_Y","Cos_Z"]
depth = pd.DataFrame(depth,columns=depthcols,dtype="float")
#print(depth)
transmitted = pd.DataFrame(transmitted,columns=transcols)
print(transmitted)

opt_data = transmitted
range_data = depth
# print(transmitted)
# h1 = transmitted.hist(column="Energy",bins=50)#
# h1 = h1[0]
# h1[0].set_xlabel("Energy(eV)")
# h1[0].set_ylabel("Counts")
# plt.show()

#Two D Histograms
# fig,ax = plt.subplots()
# cmap = "viridis"
# heatmap,xedges,yedges = np.histogram2d(transmitted["Lateral_Y"],transmitted["Lateral_Z"],bins=50)
# rmsx = np.sqrt(np.sum(np.power(transmitted["Lateral_Y"],2) )/len(transmitted["Lateral_Y"]) )
# rmsy = np.sqrt(np.sum(np.power(transmitted["Lateral_Z"],2) )/len(transmitted["Lateral_Z"]) )
# ax.set_xlabel("x ($\AA$)")
# ax.set_ylabel("y ($\AA$)")
# extent = [xedges[0], xedges[-1], yedges[0], yedges[-1]]
# minx, maxx, miny, maxy, endlim = findmax(transmitted["Lateral_Y"], transmitted["Lateral_Z"])
# ax.set_xlim((-endlim, endlim))
# ax.set_ylim((-endlim, endlim))
# tt = ax.imshow(heatmap.T, extent=extent, origin='lower', cmap=cmap, aspect="auto")
# #tt = ax.imshow(heatmap.T, extent=extent, origin='lower', cmap=cmap, aspect="equal")
# ax.set_facecolor('#440154')
# str = "$\sigma_{RMSX} =$ %.2f ($\AA$) \n$\sigma_{RMSY} =$ %.2f ($\AA$)" % (rmsx, rmsy)
# ax.plot([], [], ' ', label=str)
# ax.legend(frameon=False, facecolor='#440154', labelcolor="w")
# cbar = plt.colorbar(tt, ax=ax)
# ax.set_title("Spatial heatmap of protons leaving SRIM simulated foil \n of %s $\AA$ thickness"% (Foil_Width))
# fig.set_tight_layout(tight=True)
# fig.show()

#Max Depth Plotting
# h2,bins = np.histogram(depth["Depth_X"],bins=50)
# h2 = h2[0]
# print(bins)
# print(h2)
#
# print(len(bins))
# print(len(h2[0]))
# fig,ax = plt.subplots()
# #ax.bar(bins[:-1],h2[0])
# ax.scatter(bins,h2[0])
# fig.show()

fig,ax = plt.subplots()
#transed= depth["Depth_X"][depth["Depth_X"]>PRAL_thickness].count()
transed = depth["Depth_X"][depth["Depth_X"]>PRAL_thickness].count()
percent_transed = (transed/total_num_sim)*100
ax.hist(depth["Depth_X"],bins=50,color="royalblue",edgecolor=(0,0,0,0.4),label="Counts")
ax.set_xlabel("Range of stopped protons ($\AA$)")
ax.set_ylabel("Counts")
ax.set_title("SRIM 100 keV protons into silicon")
ax.axvline(PRAL_thickness,lw=2,color="Red",label="Suggested $\mathbf{PRAL}$ thickness for\n5 keV transmission (%i $\AA$)" % (PRAL_thickness))
ax.legend(frameon=False,title=("%.2f %% of protons travel\npast suggested thickness"%(percent_transed) ))
ax.set_xlim([0,max(depth["Depth_X"])*1.05])
fig.show()

fig,ax = plt.subplots()
transed= depth["Depth_X"][depth["Depth_X"]>Foil_Width].count()
percent_transed = (transed/total_num_sim)*100
ax.hist(depth["Depth_X"],bins=50,color="royalblue",edgecolor=(0,0,0,0.4),label="Counts")
ax.set_xlabel("Range of stopped protons ($\AA$)")
ax.set_ylabel("Counts")
ax.set_title("SRIM 100 keV protons into silicon")
ax.axvline(Foil_Width,lw=2,color="Red",label="Calculated optimal thickness for\n5 keV transmission (%i $\AA$)" % (Foil_Width))
ax.axvline(PRAL_thickness,lw=2,color="Grey",ls="--",label="$\mathbf{PRAL}$ (%i $\AA$)" % (PRAL_thickness))
ax.legend(frameon=False,title=("%.2f %% of protons travel\npast optimised thickness"%(percent_transed) ))
ax.set_xlim([0,max(depth["Depth_X"])*1.05])
fig.show()


fig,ax = plt.subplots()
trappable = transmitted["Energy"][transmitted["Energy"]<5000]#.count()
not_trappable = transmitted["Energy"][transmitted["Energy"]>5000]

print(trappable)
trappable = trappable.count()
percentage = (trappable/len(transmitted["Energy"]))*100
ax.hist(transmitted["Energy"],bins=50,color="royalblue",edgecolor=(0,0,0,0.4),label="Counts")
ax.set_xlabel("Energy of transmitted protons (eV)")
ax.set_ylabel("Counts")
ax.set_title("SRIM 100 keV protons into 8571 $\AA$ silicon")
#ax.set_title("SRIM 100 keV protons into 7945 $\AA$ silicon")
ax.axvline(5000,lw=2,color="Red",label="Maximum trappable energy")
total_trap = (total_num_sim*(percent_transed/100)*(percentage/100))
plt.rcParams['legend.title_fontsize'] = 'small'
ax.legend(frameon=False,title=("%.2f %% of transmitted protons trappable \n %2.f protons in total from 30,000 simulated"%(percentage,total_trap)))
#ax.margins(x=[0,0.05])
ax.set_xlim([0,max(transmitted["Energy"])*1.05])
fig.show()

print(percent_transed)
print(percentage)


print(len(depth["Depth_X"]))
print(trappable)
print(len(transmitted["Energy"]))
print(not_trappable.count())

#fig,ax = plt.subplot(projection="polar")
# Bootstrap resampling for errors
capture_energy = 5000
optimised_thick = 7945
n_replicas = 100
bsfrac = 0.1


print(opt_data)
earray = []
frac_array = []
per_trappable = []

depthcols = ["Ion_Number","Depth_X","Lateral_Y","Lateral_Z"]
transcols = ["Ion Numb","Atom Numb","Energy","Depth_X","Lateral_Y","Lateral_Z","Cos_X","Cos_Y","Cos_Z"]
for i in range(n_replicas):
    r = opt_data.sample(frac=bsfrac, replace=True)  # Sample the data
    x = r[r["Depth_X"] > optimised_thick]  # Taking only those that pass the thickness

    pp = range_data.sample(frac=bsfrac, replace=True)
    Nparttrans = pp["Depth_X"][pp["Depth_X"] > optimised_thick].count()
    fractrans = Nparttrans / (pp["Depth_X"].count())
    count = x["Energy"][
        x["Energy"] < capture_energy].count()  # Count those that have an energy less than the capture energy

    # tot = r["energy"].count()
    tot = x["Energy"].count()
    per_trap = count / tot * 100
    per_trappable.append(per_trap)
    per = count / tot * 100 * fractrans
    earray.append(per)
    frac_array.append(fractrans * 100)
earray = np.array(earray)
sd = earray.std()
mean = (earray.mean())
per_trappable = np.array(per_trappable)

print(sd)
N = len(earray)
# err = sd/np.sqrt(mean)
err1 = sd / np.sqrt(N)
frac_array = np.array(frac_array)
err2 = frac_array.mean() / np.sqrt(N)
err = ((err1 / mean) ** 2 + (err2 / frac_array.mean()) ** 2) ** (1 / 2)
print(err)

print(r"Average  %.2f $\pm$ %.2f" % (mean, err))
print(per_trappable.std())
print(per_trappable)
print(r"Percent trappable %.2f \pm %.2f" % (per_trappable.mean(), (per_trappable.std() / np.sqrt(N))))
print(frac_array)

print(r"frac trann trappable %.2f \pm %.2f" % (frac_array.mean(), (frac_array.std() / np.sqrt(N))))
print(frac_array)
