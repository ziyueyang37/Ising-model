import numpy as np 
import matplotlib.pyplot as plot

N = 400 #no. of spins
x = np.linspace(0.1, 4, 64)
Tc = 2/np.log(1+np.sqrt(2))  #which is nearly equals to 2.269

Mean_Collect = []
Std_Collect = []
for temp in x:

	with open("Data/Mag/Temp_%s.csv" %temp, "r") as infile:
		Data = infile.readlines() 

		Abs_Data = np.array([abs(float(x)) for x in Data[2:]])

		MPS = np.mean(Abs_Data/N)  # <|Magnetization|> per spin 
		STDPS = np.std(Abs_Data[2:]/N) # standard deviation per spin
		Mean_Collect.append(MPS)
		Std_Collect.append(STDPS)

#Mean field
MF = []
s = 1
for temp in x:
    for _ in range(1000):
        s = np.tanh(Tc/temp * s)

    MF.append(s)


#Coding for Osager Solution's for Finding Curie Temperature.


OM = [] #Osager's Magnetization
for temp in x:
	if float(temp) < Tc:
		O =(1-(np.sinh(np.log(1+np.sqrt(2))*(Tc/temp)))**(-4))**(1/8)
		OM.append(O)
	elif float(temp) > Tc:
		O = 0
		OM.append(O)

f1 = plot.figure()
plot.plot(x/Tc,MF, color = "darkblue", marker = "*", label = "Mean Field Approximation")
plot.grid()
plot.title("Absolute Magnetization Mean Field Approximation")
plot.legend(loc = "best")
plot.xlabel("Reduced Temperature T/Tc")
plot.ylabel("Absolute Magnetization per spin")
plot.savefig('./figures/Magplot_mf.png')

f2 = plot.figure()
plot.plot(x/Tc, OM, color = "salmon", marker = "*", label = "Onsager's formula")
plot.grid()
plot.title("Absolute Magnetization Onsager's formula")
plot.legend(loc = "best")
plot.xlabel("Reduced Temperature T/Tc")
plot.ylabel("Absolute Magnetization per spin")
plot.savefig('./figures/Magplot_om.png')

f3 = plot.figure()
plot.errorbar(x/Tc,Mean_Collect, yerr = Std_Collect, color= "forestgreen", fmt = "o", label = "Simulation Data")
plot.grid()
plot.title("Absolute Magnetization Simulation")
plot.legend(loc = "best")
plot.xlabel("Reduced Temperature T/Tc")
plot.ylabel("Absolute Magnetization per spin")
plot.savefig('./figures/Magplot_sim.png')


