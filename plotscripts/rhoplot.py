from yt.mods import *
from math    import sqrt
from os      import system

# Defining new field to plot
def _velt(field,data):
	vx2 = (data["velx"])**2
	vy2 = (data["vely"])**2
	vz2 = (data["velz"])**2
	vt  = (vx2 + vy2 + vz2)**0.5
	return (vt)
add_field("velt",function=_velt,units=r"\rm{cm}/\rm{s}")

cmd = 'mkdir -p frames'
os.system(cmd)

fn     = "wdwd_hdf5_plt_cnt_"
x      = 0.0e+10
y      = 0.0e+10
z      = 0.0e+10

L = ["velt"]
# Loop over all plot variables
for variable in L:
  os.system(cmd+"/"+variable)
  j = 0
  # Loop over desired files
  for i in range(0,10):
    if -1  < i < 10  : ni = "000" + str(i)
    if  9  < i < 100 : ni = "00"  + str(i)
    if 100 < i < 1000: ni = "0"   + str(i)
    dmp = " >> frames/"+variable+"/log"
    na = fn+ni
    pf = load(na)
    dd = pf.h.all_data()
    # Find extrema
    ext   = dd.quantities["Extrema"](variable)[0]
    min   = ext[0]
    max   = ext[1]
    ext   = dd.quantities["Extrema"]("velx")[0]
    vxmin = ext[0]
    vxmax = ext[1]
    ext   = dd.quantities["Extrema"]("vely")[0]
    vymin = ext[0]
    vymax = ext[1]
    ext   = dd.quantities["Extrema"]("velz")[0]
    vzmin = ext[0]
    vzmax = ext[1]
    ext   = dd.quantities["Extrema"]("velt")[0]
    vtmin = ext[0]
    vtmax = ext[1] 
    # Echo them back
    os.system("echo time ="+str(ni)+dmp)
    os.system("echo vamin="+str(min)+dmp)
    os.system("echo vamax="+str(max)+dmp) 
    os.system("echo vxmin="+str(vxmin)+dmp)
    os.system("echo vxmax="+str(vxmax)+dmp) 
    os.system("echo vymin="+str(vymin)+dmp)
    os.system("echo vymax="+str(vymax)+dmp) 
    os.system("echo vzmin="+str(vzmin)+dmp)
    os.system("echo vzmax="+str(vzmax)+dmp)  
    os.system("echo vtmin="+str(vtmin)+dmp)
    os.system("echo vtmax="+str(vtmax)+dmp)  
    os.system("echo       "+dmp)   
    # Let's project
    center = [x,y,z]
    pc = PlotCollection(pf,center)
    p2 = pc.add_slice(variable,"z")
    p2.modify["velocity"](normalize=True)
    pc.save("frames/"+variable+"/"+na)
