# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.10.0
#   kernelspec:
#     display_name: Python 3
#     language: python
#     name: python3
# ---

# +
###################################################################################################
# Read ASCII text output file from RAMS SCM into numpy arrays and plot some vertical profiles.
###################################################################################################
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
# Setting font parameters for title and axes
font = {'family' : 'verdana',
        'size'   : 16,
        'weight' : 'bold'}
mpl.rc('font', **font)
mpl.rcParams['axes.linewidth'] = 2  #Set the value globally
# Set some fonts here for easy changing
font_axis      = 22
font_suptitle  = 18
font_xylabel   = 22
font_gridlabel = 11
font_colorbar  = 12
font_legend    = 12
font_title     = 18
###################################################################################################
# Function to stop at a certain point in the jupyter notebook cell.
# Use something like 'raise sms.StopExecution' to stop code in a cell when using import
# method above.
class StopExecution(Exception):
    def _render_traceback_(self):
        pass
###################################################################################################
# Set whether to save or show the figure image
plot_to_file=True

# Set variable to plot and in/out directories
varname='rrp'
infile='scm.sedimentation-test.in'
outfile='scm.sedimentation-test.out'

# Read the vertical scalar levels file 
f=open(f'{infile}/zt.txt','r')
line = f.readlines()
# Get number of vertical levels (top line is a header to subtract)
numlevs = len(line)-1
# Create array for the levels
ztlev = np.zeros((numlevs),dtype='float32')
# Read the levels to the array
for i in range(len(line)-1):
    value = line[i+1].split()
    ztlev[i] = float(value[-1])
f.close()
ztlev = ztlev / 1000.
#print('Number of vertical levels:',numlevs,ztlev)

# First get the array of times from an output file
f=open(f'{outfile}/{varname}.txt','r')
numtimes=0
while True:
    line = f.readline()
    if not line:
        break
    value = line.split()
    if value[0] == 'Elapsed':
        numtimes+=1
        if numtimes == 1:
            times = np.array([float(value[-1])])
        else:
            times = np.append(times,float(value[-1]))
f.close()
#print('Number of times:',numtimes,times)

# Get the actually data values and place in array that we can plot
vari = np.zeros((numtimes,numlevs),dtype='float32')
tim=-1
lev=-1
count=0
f=open(f'{outfile}/{varname}.txt','r')
while True:
    line = f.readline()
    if not line:
        break
    count+=1
    value = line.split()
    if value[0] == 'Elapsed':
        tim += 1
        lev = -1
    else:
        lev += 1
        vari[tim,lev]=value[-1]
        #print(tim,lev,times[tim],ztlev[lev],vari[tim,lev])
f.close()

if varname == 'rrp':
    vari = vari * 1000. #convert kg/kg to g/kg
    xlabel = 'Rain Mixing Ratio ($g/kg$)'

###################################################################################################
for tim in range(numtimes):
    fig,ax=plt.subplots(nrows=1,ncols=1,figsize=(7,7))
    ax.set_ylabel('Altitude (km)',fontsize=font_xylabel, fontweight='bold')
    ax.set_xlabel(xlabel,fontsize=font_xylabel,fontweight='bold')
    plt.xticks(fontsize=16)
    plt.yticks(fontsize=16)
    ax.set_ylim(ymin=0.0)
    ax.set_yticks(np.arange(0, 5, 1))
    xlimmax = np.max(vari)
    xlimmin = np.min(vari)
    ax.set_xlim([xlimmin,xlimmax])
    ax.set_title(f'{times[tim]} Sec Elapsed Time',fontsize=font_title,fontweight='bold')
    ax.grid()
    ax.plot(vari[tim],ztlev,color='blue',linestyle='solid',linewidth=3)
    fig.patch.set_facecolor('white')
    fig.tight_layout()
    plotname=f'./scm.sedimentation-test.out/Plt.VP.{varname}.t{tim:03d}.png'
    if plot_to_file is True:
        print(tim,plotname)
        fig.savefig(plotname,bbox_inches='tight',pad_inches=0.2,dpi=200)
        plt.close(fig)
    else:
        plt.show()

print('FINISHED')
###################################################################################################
# -
