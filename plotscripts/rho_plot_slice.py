from yt.mods import *
from os      import system

# Commenting this for the sake of readability
# Hi, it's ok if you don't know the full syntax for this
# I googled how to do a lot of things in python

# This script creates slices of your simulation 
# along the z axis as a normal vector

# fn is the prefix of tjhe input files
# fo is the prefix of your output plots
# varlist is the list of variables you want to plot
# tmin is the minimum time inex
# tmax is the maximum time index
fn      = "wdwd_hdf5_plt_cnt_"
fo      = "wdwd_plt_"
varlist = ["velt"]
tmin = 0
tmax = 25

# x,y are the center of your simulation
# zcen is the central slice, kmax the maximum distance from zcen you want to slice to
# kmax is the number of slices between the zcen slice and the zmax slice
# slices will also be created between zcen and -zmax as well
x    = 0.0e+10
y    = 0.0e+10
zcen = 0.0e+10
zmax = 2.461e9
kmax = 8
dz   = (zmax-zcen) / float(kmax)

# Defining the norm of the velocity and letting YT know about it
def _velt(field,data):
	vx2 = (data["velx"])**2
	vy2 = (data["vely"])**2
	vz2 = (data["velz"])**2
	vt  = (vx2 + vy2 + vz2)**0.5
	return (vt)
add_field("velt",function=_velt,units=r"\rm{cm}/\rm{s}")

# Create working directory
cmd = 'mkdir -p frames'
os.system(cmd)

# Loop over all plot variables
tdum = tmax
switch1 = 0
switch2 = 0
for i in range(tmin,tdum+1):
  if -1  < i < 10  : ni = "000" + str(i)
  if  9  < i < 100 : ni = "00"  + str(i)
  if 100 < i < 1000: ni = "0"   + str(i)\
  na = fn+ni
  no = fo+ni
  pf = load(na)
  dd = pf.h.all_data()
  if switch2 == 1: switch2 = 1
  for variable in varlist:
    # Opening a log file and creating a directory if needed
    os.system(cmd+"/"+variable)
    if switch1 == 0: 
      lg = open("frames/"+variable+"/slices_log","w")
      lg.write("# tindex "+variable+"_max maxloc\n")
    lg = open("frames/"+variable+"/slices_log","a")
    # Find extrema
    ext   = dd.quantities["Extrema"](variable)[0]
    min   = ext[0]
    max   = ext[1]
    thing = pf.h.find_max(variable)
    dum11 = thing[1]
    xxmax = dum11[0]
    yymax = dum11[1]
    zzmax = dum11[2]
    ext   = dd.quantities["Extrema"]("velt")[0]
    vtmin = ext[0]
    vtmax = ext[1] 
    # Logging
    lg.write(ni+"   "+str(max)+"   "+str(xxmax)+"   "+str(yymax)+"   "+str(zzmax)+"\n")
    # Let's slice this thing up 
    os.system(cmd+"/"+variable+"/"+"slices_"+str(ni))
    for k in range (-kmax,kmax+1): 
      if -1  < k+kmax < 10  : nk = "000" + str(k+kmax)
      if  9  < k+kmax < 100 : nk = "00"  + str(k+kmax)
      if 100 < k+kmax < 1000: nk = "0"   + str(k+kmax)
      z = zcen + float(k)*dz
      center = [x,y,z]
      pc = PlotCollection(pf,center)
      p2 = pc.add_slice(variable,"z")
      p2.modify["velocity"](normalize=True)
      pc.save("frames/"+variable+"/slices_"+str(ni)+"/"+no+"_"+nk)
  switch2 = 1
