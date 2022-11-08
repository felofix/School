import numpy as np
import matplotlib.pyplot as plt

navn_tekst = np.array(['vectors/eigenvector10.txt', 'vectors/eigenvector11.txt','vectors/eigenvector12.txt'])
navn_tekst_anal = np.array(['vectors/eigenvectoranal10.txt', 'vectors/eigenvectoranal11.txt','vectors/eigenvectoranal12.txt'])

navn_tittel = np.array(["1. smallest", "2. smallest", "3. smallest"])
linestyles = np.array(['dotted', 'dashed', 'dashdot'])
autumn_colors = np.array(['#5EFF36', '#1C24FF', '#FF2B47'])

for i in range(len(navn_tekst)):
    data = np.loadtxt(navn_tekst[i])
    data_anal = np.loadtxt(navn_tekst_anal[i])
    x = data[:,0]
    xa = data[:,0]
    u_x = data[:,1]
    u_ax = data[:,1]
    u_x[0] = 0
    u_x[-1] = 0
    u_ax[0] = 0
    u_ax[-1] = 0
    plt.plot(x, u_x, color = autumn_colors[i], label = f"Approx {navn_tittel[i]}")
    plt.plot(xa, u_ax, color = 'grey', linestyle = linestyles[i], label = f"Analytical {navn_tittel[i]}")
   
plt.title("Three lowest eigenvectors for buckling beam problem with n = 100")
plt.xlabel("$\hat{x}$")
plt.ylabel("$ u(\hat{x}$)")
plt.legend(prop={'size': 8})
plt.savefig("Plot for the three smallest eigenvalues_n100.pdf")

    
