import numpy as np 
import matplotlib.pyplot as plot

N = 400/10000 #no. of spins
x = np.linspace(0.1, 4, 64)
dT = (4-0.1)/64
Tc = 2/np.log(1+np.sqrt(2)) 

Mean_Collect = []
Std_Collect = []
for temp in x:

	with open("Data/cv/Temp_%s.csv" %temp, "r") as infile:
		Data = infile.readlines() 

		Abs_Data = np.array([abs(float(x)) for x in Data[2:]])

		MPS = np.mean(Abs_Data/N/10)  # <|Magnetization|> per spin 
		STDPS = np.std(Abs_Data[2:]/N/10) # standard deviation per spin
		Mean_Collect.append(MPS)
		Std_Collect.append(STDPS)
MFCv=[]
MFE = []
s = 1
for temp in x:
    for _ in range(1000):
        s = np.tanh(Tc/temp * s)
    MFE.append(Tc*s**2)
for i in range(1,len(MFE)):
    MFCv.append(-(MFE[i]-MFE[i-1])/dT**2/20)
MFCv.append(0)


OMCv=[]
OME = []
OM=[]
for temp in x:
    O = 0
    if float(temp) < Tc:
        O =(1-(np.sinh(np.log(1+np.sqrt(2))*(Tc/temp)))**(-4))**(1/8)
    OME.append(Tc*O**2)



for i in range(1,len(OME)):
    OMCv.append(-(OME[i]-OME[i-1])/dT**2/80)
OMCv.append(0)


plot.errorbar(x/Tc,Mean_Collect, yerr = Std_Collect, color= "maroon", fmt = "o")
plot.plot(x/Tc, Mean_Collect, marker = "*", color="maroon", label="Simulation")
plot.title("Specific Heat Simulation Result")
plot.legend(loc = "best")
plot.grid()
plot.xlabel("Reduced Temperature Tc/T")
plot.ylabel("Specific Heat")
plot.savefig('./figures/cvplot_sim.png')

f2 = plot.figure()
plot.plot(x/Tc, MFCv, marker="*", color="steelblue", label="Mean Field Appro.")
plot.title("Specific Heat Mean Field Approximation")
plot.legend(loc = "best")
plot.grid()
plot.xlabel("Reduced Temperature Tc/T")
plot.ylabel("Specific Heat")
plot.savefig('./figures/cvplot_mf.png')

f3 = plot.figure()
plot.plot(x/Tc, OMCv, marker="*", color="slategray", label="Onsager's formula")
plot.title("Specific Heat Onsager's formula")
plot.legend(loc = "best")
plot.grid()
plot.xlabel("Temperature")
plot.ylabel("Specific Heat")
plot.savefig('./figures/cvplot_om.png')

