from yt.mods import *
import os
cmd = 'mkdir -p frames'
os.system(cmd)

fn     = "wdwd_hdf5_plt_cnt_"
x      = 0.0e+10
y      = 0.0e+10
z      = 0.0e+10
center = [x,y,z]

L = ["dens"] #["dens","gpot","lapl","gpsi","Grav_Potential"]
M = [0,1,2,3]
for variable in L:
  os.system(cmd+"/"+variable)
  for i in M:
    if -1  < i < 10  : ni = "000" + str(i)
    if  9  < i < 100 : ni = "00"  + str(i)
    if 100 < i < 1000: ni = "0"   + str(i)
    na = fn+ni
    pf = load(na)
    dd = pf.h.all_data()
    mi,ma = dd.quantities["Extrema"](variable)[0]
    pc = PlotCollection(pf,center)
    p2 = pc.add_slice(variable,"z")
    pc.set_zlim(mi,ma)
    pc.save("frames/"+variable+"/"+na)
