import numpy as np
import matplotlib.pyplot as plt


SEX_full=np.genfromtxt('o.YPP-X_probe_order_2',comments="#")

fig=plt.figure(figsize=(10,10))

#pmVm1toCGS=2.38721e-09
pmVm1toCGS=1
plt.title('SHG in hBN: SEX')
plt.xlabel('eV')
plt.ylabel('$|\chi^2_{zxy}|$ (pm/V)',fontsize=18)
#plt.xlim(1.,5)
#plt.ylim(0,1400)
plt.plot(SEX_full[:,0], np.sqrt(SEX_full[:,1]**2+SEX_full[:,2]**2)/pmVm1toCGS,color="grey", linestyle="--",label="$|\chi^2|$ (20 freqs ABINIT)")

plt.legend(fontsize=18)
plt.savefig("hBN-SHG.png",format="png")
plt.show()

