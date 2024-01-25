import numpy as np
D=1
Length=np.linspace(0,10,101)
         
                
#Create stenosis tags for each of the centerline point
                
Loc_stenosis=Length[-1]/100.*50 #location of stenosis along CL

s0=0.25; x_=0.0;L_stenosis=2*D #Varghese et al. stenosis parameters
                
for i in range(0,100):
	if (Length[i]>=(Loc_stenosis-D)) and (Length[i]<=(Loc_stenosis+D)):
		x_+=Length[i]-Length[i-1]
		x0=Loc_stenosis 
		Dia=0.5*D*(1-s0*(1+np.cos(2*np.pi*(x_-x0)/L_stenosis)))
		print (Length[i],Dia)
