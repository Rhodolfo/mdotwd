from yt.mods import *
import os
cmd = 'mkdir -p frames'
os.system(cmd)

fn     = "rho_hdf5_plt_cnt_"
x      = 0.0e+10
y      = 0.0e+10
z      = 0.0e+10
center = [x,y,z]

L = ["dens"]
for variable in L:
  os.system(cmd+"/"+variable)
  j = 0
  for i in range(350,444):
    if -1  < i < 10  : ni = "000" + str(i)
    if  9  < i < 100 : ni = "00"  + str(i)
    if 100 < i < 1000: ni = "0"   + str(i)
    na = fn+ni
    pf = load(na)
    dd = pf.h.all_data()
    if j == 0:
      mi, ma = dd.quantities["Extrema"](variable)[0]
      if variable == "dens":
        mi = 1.0e-15
      j = 1
    pc = PlotCollection(pf,center)
    p2 = pc.add_projection(variable,"z")
    pc.set_zlim(mi,ma)
    pc.save("frames/"+variable+"/"+na)
