# run this as python nanotest.py > nanotest.log &

import os
import nanogen as ng
import myutil as my

pD = my.setParamDefaults()
pD['Duration'] = 10000
pD['Folder'] = 'Nanotube With CUDA Test'
pD['m'] = 5
pD['n'] = 5

force = [ 0.01, 0.02, 0.04, 0.08, 0.16 ]
temperature = [ 80 ]
for f in force:
    for t in temperature:
        pD['Force (pN)'] = f
        pD['Temperature (K)'] = t
        
        for r in range(1, 2):
            pD['File Name'] = 'Run-' + str(r)
            ng.run(pD)