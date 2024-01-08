prepgrid.exe
p.exe
copy param_br param_br.dat
copy grid_new.dat grid_bnd.dat
copy nodal_g nodal_g.dat
erase nodal_g
erase grid_new.dat
erase param_br
prepdiff.exe
copy paratemp.dat param_g.dat
erase paratemp.dat
prepsurf.exe
stomata.exe
