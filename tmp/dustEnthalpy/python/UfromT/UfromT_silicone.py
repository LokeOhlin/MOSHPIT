import numpy as np
import matplotlib.pyplot as plt


def UfromT(T):
    if (T < 50):
        U = 1.4*1e3/3 * T**3
    else:
        U = 1.4*1e3/3* (50)**3
        if (T < 150): 
            U += 2.2*1e4 / (2.3) * ( T**(2.3) - 50**(2.3))
        else: 
            U += 2.2*1e4 / (2.3) * ( 150**(2.3) - 50**(2.3))
            if (T < 500):
                U += 4.8*1e5/(1.68) * (T**(1.68) - 150**(1.68))
            else:
                U += 4.8*1e5/(1.68) * (500**(1.68) - 150**(1.68))
                U += 3.41*1e7*(T - 500)
    return U 

def TfromU(U):
    Ulim1 = UfromT(50)
    Ulim2 = UfromT(150)
    Ulim3 = UfromT(500)
    print(Ulim1)
    print(Ulim2)
    print(Ulim3)
    Uint = 0
    if U < Ulim1:
        return (U/(1.4e3 / 3))**(1/3)
    else:
        Uint = 1.4*1e3/3* (50)**3
        print("Uint1",1.4*1e3/3* (50)**3)

        if(U < Ulim2):
            return ((U-Uint)/(2.2*1e4 / (2.3)) + 50**2.3)**(1/2.3)
        else: 
            Uint = Uint + 2.2*1e4 / (2.3) * ( 150**(2.3) - 50**(2.3))
            print("Uint2",Uint)
            if(U < Ulim3):
                return ((U-Uint)/(4.8*1e5/(1.68)) + 150**1.68)**(1/1.68)
            else:
                Uint = Uint + 4.8*1e5/(1.68) * (500**(1.68) - 150**(1.68))
                print("Uint3",Uint)
                return ((U-Uint)/(3.41*1e7) + 500)


temps = np.logspace(0, 4)
temps_inv = np.zeros(temps.shape)
maxDiff = 0
for i in range(len(temps)):
    U = UfromT(temps[i])
    temps_inv[i] = TfromU(U)
    maxDiff = max(maxDiff, abs(temps_inv[i] - temps[i]))
print(maxDiff)
plt.plot(temps, temps_inv)
plt.xscale('log')
plt.yscale('log')
plt.show()
