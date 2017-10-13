# -*- coding: utf-8 -*-
"""
Created on Fri Dec  9 20:15:03 2016

@author: Robert
"""

import numpy as np
import matplotlib.pyplot as plt
from random import random
from math import log
from copy import copy
import seaborn as sns

import pandas


class Error(Exception):
   """Base class for other exceptions"""
   pass

class NegativeRate(Error):
   """Negative Birth or Death Rate"""
   pass


data=[]

replicates = 200

cutoff=1000

#varFactor is birth + death for all cell types in all drugs:

varFactor = 0.1

g_x_odd = -0.05 #Initial treatment (eg ALK TKI)
g_y_odd = 0.005

b_x_odd = (g_x_odd+varFactor)/2
b_y_odd = (g_y_odd+varFactor)/2

d_x_odd = (-g_x_odd+varFactor)/2
d_y_odd = (-g_y_odd+varFactor)/2

g_x_even = 0.005 
g_y_even = -0.05 

b_x_even = (g_x_even+varFactor)/2
b_y_even = (g_y_even+varFactor)/2

           
d_x_even = (-g_x_even+varFactor)/2
d_y_even = (-g_y_even+varFactor)/2

trans1 = 0.0001 
trans2 = trans1

relativePopInit = [9000,1000];

#Check for negative rates:
for parameter in [b_x_odd,b_y_odd,d_x_odd,d_y_odd,b_y_even,b_x_even,d_x_even,d_y_even]:
    if parameter < 0:
        raise NegativeRate

ApB=relativePopInit[1]/relativePopInit[0]
        
T_initial = np.log(((g_y_odd-g_x_odd)*(g_x_even-g_x_odd)+trans1*(g_y_odd+g_x_even-g_x_odd-g_y_even))/((trans1+ApB*(trans1+g_y_odd-g_x_odd))*(g_y_odd-g_y_even)))/(trans1+g_y_odd-g_x_odd)


k=(((g_y_odd-g_x_odd)*(g_x_even-g_x_odd)+trans1*(g_y_odd+g_x_even-g_x_odd-g_y_even))*(g_y_odd-g_y_even))/(((g_y_odd-g_y_even)*(g_x_even-g_y_even)+trans1*(g_y_odd+g_x_even-g_x_odd-g_y_even))*(g_x_even-g_x_odd))

T_treat_1 = 10.0
T_treat_2 = T_treat_1/k

T_treat = T_treat_1+T_treat_2

#T_initial -> T_treat_1 -> T_treat_2 -> T_treat_1

print(T_initial,T_treat_1,T_treat_2)



    





relativePop = relativePopInit


          
#-----------------Stoch----------------------------------------

def update(state):
    
    rt,re = random(), random()
    
    timeAfter = (state[2]-T_initial)/T_treat-int((state[2]-T_initial)/T_treat)
    
    
    if timeAfter>(T_treat_1/T_treat) or (state[2]<T_initial):
        a =(b_x_odd+d_x_odd+trans1)*state[0] + (b_y_odd+d_y_odd)*state[1]  
        
        if b_x_odd*state[0]/a>re:
            state[0] += 1
        elif (d_x_odd*state[0]+b_x_odd*state[0])/a>re:
            state[0] -= 1
        elif (b_y_odd*state[1]+d_x_odd*state[0]+b_x_odd*state[0])/a>re:
            state[1] += 1
        elif (d_y_odd*state[1]+b_y_odd*state[1]+d_x_odd*state[0]+b_x_odd*state[0])/a>re:
            state[1] -= 1
        else:
            state[0] -= 1
            state[1] += 1
    else:
        a =(b_x_even+d_x_even)*state[0] + (b_y_even+d_y_even+trans2)*state[1]
        
        if b_x_even*state[0]/a>re:
            state[0] += 1
        elif (d_x_even*state[0]+b_x_even*state[0])/a>re:
            state[0] -= 1
        elif (b_y_even*state[1]+d_x_even*state[0]+b_x_even*state[0])/a>re:
            state[1] += 1
        elif (d_y_even*state[1]+b_y_even*state[1]+d_x_even*state[0]+b_x_even*state[0])/a>re:
            state[1] -= 1
        else:
            state[0] += 1
            state[1] -= 1
             
    state[2] = state[2] - (1/a)*log(rt)
        
    
    return state



replicate = 0

extinctionTimes=[]

data=[]



while (replicate < replicates):   
    
    extinct=False;
    
    state = [relativePopInit[0], relativePopInit[1], 0.0] #initial conditions for state vector(0=sensitive, 1=resistant, 2=time)    
    
    life_history = []
    life_history.append(state)
    
    while (extinct==False) & (state[2]<cutoff):
    
        #print "Time: ", state[2]," Sensitive: ", state[0]," Resistant: ", state[1]
            
        new_state = copy(state)

        state = update(new_state)    
        
        life_history.append(state)
        
        if ((state[0]==0) & (state[1]==0)):
            extinct=True
            extinctionTimes.append(state[2])
            
    
    life_history.append([0,0,cutoff])
    
    data.append(life_history)
    
    replicate += 1
    
#-----------------ODE----------------------------------------


data0=[];
data1=[];

timeStep=0.01;

for k in range(0, int(cutoff/timeStep)):
    
    if (float(k)*timeStep).is_integer():
        data0.append(relativePop[0]);
        data1.append(relativePop[1]);
        
    r01=trans1;
    r10=trans2;    
    
    
    timeAfter = (k*timeStep-T_initial)/T_treat-int((k*timeStep-T_initial)/T_treat)
    
    if timeAfter>(T_treat_1/T_treat) or (k*timeStep<T_initial):
        g0=g_x_odd;
        g1=g_y_odd;
        r01=trans1;
        r10=0.0;
    else:
        g0=g_x_even;
        g1=g_y_even;
        r01=0.0;
        r10=trans2;
    
    
        
    relativePop[0]=relativePop[0]+(relativePop[0]*g0+relativePop[1]*r10-relativePop[0]*r01)*timeStep;
    relativePop[1]=relativePop[1]+(relativePop[1]*g1+relativePop[0]*r01-relativePop[1]*r10)*timeStep;

totalODE=[]

for i,na in enumerate(data0):
    totalODE.append(data0[i]+data1[i])    
    

                   
#------------------Plots--------------------------

preDf=[]


#Bin the data into replicates:

for replicateNumber,replicateData in enumerate(data):
    
    binData=[[[] for _ in range(cutoff+1)],[[] for _ in range(cutoff+1)]]
    
    for i,repData in enumerate(replicateData):
        
        binData[0][int(replicateData[i][2])].append(int(replicateData[i][0]))
        binData[1][int(replicateData[i][2])].append(int(replicateData[i][1]))
        
    for i,binList in enumerate(binData[0]):
        
        tempData=[0,0,0]        
        
        if len(binData[0][i])==0:
            tempData[0]=binData[0][i-1][len(binData[0][i-1])-1]
            binData[0][i].append(binData[0][i-1][len(binData[0][i-1])-1])
        else:
            tempData[0]=np.mean(binData[0][i])
        if len(binData[1][i])==0:    
            tempData[1]=binData[1][i-1][len(binData[1][i-1])-1]
            binData[1][i].append(binData[1][i-1][len(binData[0][i-1])-1])
        else:
            tempData[1]=np.mean(binData[1][i])
        
        tempData[2]=sum(tempData)
        
        preDf.append([float(i),'$A_R$',str(replicateNumber),tempData[1]]) #change order here to change order of cells on tsplot
        preDf.append([float(i),'$B_R$',str(replicateNumber),tempData[0]])
        
        preDf.append([float(i),'Total',str(replicateNumber),tempData[2]])

odepreDf=[]
odepreDf2=[]

for k in range(cutoff):
    odepreDf.append([float(k),'ODE',replicates,data0[k]])
    
for k in range(cutoff):
    odepreDf2.append([float(k),'ODE2',replicates,data1[k]]) 
    
for k in range(cutoff):
    odepreDf2.append([float(k),'ODE3',replicates+1,totalODE[k]]) 

            
df=pandas.DataFrame(preDf,columns=['Time (days)', 'Cell Type','subject', 'Cell Number'])

odeDf=pandas.DataFrame(odepreDf,columns=['Time (days)', 'Cell Type','subject', 'Cell Number'])
odeDf2=pandas.DataFrame(odepreDf2,columns=['Time (days)', 'Cell Type','subject', 'Cell Number'])

#tsplot:

#sns.set_style("whitegrid", {'axes.grid' : False})

sns.set(font_scale=1.5)

ax = sns.tsplot(time="Time (days)", value="Cell Number",unit="subject", condition="Cell Type",data=df,err_style="unit_traces")

ax = sns.tsplot(time="Time (days)", value="Cell Number",unit="subject", condition="Cell Type",data=odeDf,color='k',linestyle='--',legend=True)
ax = sns.tsplot(time="Time (days)", value="Cell Number",unit="subject", condition="Cell Type",data=odeDf2,color='k',linestyle='--',legend=False)

#ax.set(yscale="log")

ax.grid(color='k',linestyle=':')

ax.set_xlim(0,200)
#ax.set_ylim(0,500)




#print(len(df))



#print(extinctionTimes)

data01=extinctionTimes

#######################################################################################################################################
#######################################################################################################################################
varFactor = 0.2

g_x_odd = -0.05 #Initial treatment (eg ALK TKI)
g_y_odd = 0.005

b_x_odd = (g_x_odd+varFactor)/2
b_y_odd = (g_y_odd+varFactor)/2

d_x_odd = (-g_x_odd+varFactor)/2
d_y_odd = (-g_y_odd+varFactor)/2

g_x_even = 0.005 
g_y_even = -0.05 

b_x_even = (g_x_even+varFactor)/2
b_y_even = (g_y_even+varFactor)/2

           
d_x_even = (-g_x_even+varFactor)/2
d_y_even = (-g_y_even+varFactor)/2

trans1 = 0.0001 
trans2 = trans1

relativePopInit = [9000,1000];

#Check for negative rates:
for parameter in [b_x_odd,b_y_odd,d_x_odd,d_y_odd,b_y_even,b_x_even,d_x_even,d_y_even]:
    if parameter < 0:
        raise NegativeRate

ApB=relativePopInit[1]/relativePopInit[0]
        
T_initial = np.log(((g_y_odd-g_x_odd)*(g_x_even-g_x_odd)+trans1*(g_y_odd+g_x_even-g_x_odd-g_y_even))/((trans1+ApB*(trans1+g_y_odd-g_x_odd))*(g_y_odd-g_y_even)))/(trans1+g_y_odd-g_x_odd)


k=(((g_y_odd-g_x_odd)*(g_x_even-g_x_odd)+trans1*(g_y_odd+g_x_even-g_x_odd-g_y_even))*(g_y_odd-g_y_even))/(((g_y_odd-g_y_even)*(g_x_even-g_y_even)+trans1*(g_y_odd+g_x_even-g_x_odd-g_y_even))*(g_x_even-g_x_odd))

T_treat_1 = 10.0
T_treat_2 = T_treat_1/k

T_treat = T_treat_1+T_treat_2

#T_initial -> T_treat_1 -> T_treat_2 -> T_treat_1

print(T_initial,T_treat_1,T_treat_2)



    





relativePop = relativePopInit


          
#-----------------Stoch----------------------------------------

def update(state):
    
    rt,re = random(), random()
    
    timeAfter = (state[2]-T_initial)/T_treat-int((state[2]-T_initial)/T_treat)
    
    
    if timeAfter>(T_treat_1/T_treat) or (state[2]<T_initial):
        a =(b_x_odd+d_x_odd+trans1)*state[0] + (b_y_odd+d_y_odd)*state[1]  
        
        if b_x_odd*state[0]/a>re:
            state[0] += 1
        elif (d_x_odd*state[0]+b_x_odd*state[0])/a>re:
            state[0] -= 1
        elif (b_y_odd*state[1]+d_x_odd*state[0]+b_x_odd*state[0])/a>re:
            state[1] += 1
        elif (d_y_odd*state[1]+b_y_odd*state[1]+d_x_odd*state[0]+b_x_odd*state[0])/a>re:
            state[1] -= 1
        else:
            state[0] -= 1
            state[1] += 1
    else:
        a =(b_x_even+d_x_even)*state[0] + (b_y_even+d_y_even+trans2)*state[1]
        
        if b_x_even*state[0]/a>re:
            state[0] += 1
        elif (d_x_even*state[0]+b_x_even*state[0])/a>re:
            state[0] -= 1
        elif (b_y_even*state[1]+d_x_even*state[0]+b_x_even*state[0])/a>re:
            state[1] += 1
        elif (d_y_even*state[1]+b_y_even*state[1]+d_x_even*state[0]+b_x_even*state[0])/a>re:
            state[1] -= 1
        else:
            state[0] += 1
            state[1] -= 1
             
    state[2] = state[2] - (1/a)*log(rt)
        
    
    return state



replicate = 0

extinctionTimes=[]

data=[]



while (replicate < replicates):   
    
    extinct=False;
    
    state = [relativePopInit[0], relativePopInit[1], 0.0] #initial conditions for state vector(0=sensitive, 1=resistant, 2=time)    
    
    life_history = []
    life_history.append(state)
    
    while (extinct==False) & (state[2]<cutoff):
    
        #print "Time: ", state[2]," Sensitive: ", state[0]," Resistant: ", state[1]
            
        new_state = copy(state)

        state = update(new_state)    
        
        life_history.append(state)
        
        if ((state[0]==0) & (state[1]==0)):
            extinct=True
            extinctionTimes.append(state[2])
            
    
    life_history.append([0,0,cutoff])
    
    data.append(life_history)
    
    replicate += 1
    
#-----------------ODE----------------------------------------


data0=[];
data1=[];

timeStep=0.01;

for k in range(0, int(cutoff/timeStep)):
    
    if (float(k)*timeStep).is_integer():
        data0.append(relativePop[0]);
        data1.append(relativePop[1]);
        
    r01=trans1;
    r10=trans2;    
    
    
    timeAfter = (k*timeStep-T_initial)/T_treat-int((k*timeStep-T_initial)/T_treat)
    
    if timeAfter>(T_treat_1/T_treat) or (k*timeStep<T_initial):
        g0=g_x_odd;
        g1=g_y_odd;
        r01=trans1;
        r10=0.0;
    else:
        g0=g_x_even;
        g1=g_y_even;
        r01=0.0;
        r10=trans2;
    
    
        
    relativePop[0]=relativePop[0]+(relativePop[0]*g0+relativePop[1]*r10-relativePop[0]*r01)*timeStep;
    relativePop[1]=relativePop[1]+(relativePop[1]*g1+relativePop[0]*r01-relativePop[1]*r10)*timeStep;

totalODE=[]

for i,na in enumerate(data0):
    totalODE.append(data0[i]+data1[i])    
    

                   
#------------------Plots--------------------------

preDf=[]


#Bin the data into replicates:

for replicateNumber,replicateData in enumerate(data):
    
    binData=[[[] for _ in range(cutoff+1)],[[] for _ in range(cutoff+1)]]
    
    for i,repData in enumerate(replicateData):
        
        binData[0][int(replicateData[i][2])].append(int(replicateData[i][0]))
        binData[1][int(replicateData[i][2])].append(int(replicateData[i][1]))
        
    for i,binList in enumerate(binData[0]):
        
        tempData=[0,0,0]        
        
        if len(binData[0][i])==0:
            tempData[0]=binData[0][i-1][len(binData[0][i-1])-1]
            binData[0][i].append(binData[0][i-1][len(binData[0][i-1])-1])
        else:
            tempData[0]=np.mean(binData[0][i])
        if len(binData[1][i])==0:    
            tempData[1]=binData[1][i-1][len(binData[1][i-1])-1]
            binData[1][i].append(binData[1][i-1][len(binData[0][i-1])-1])
        else:
            tempData[1]=np.mean(binData[1][i])
        
        tempData[2]=sum(tempData)
        
        preDf.append([float(i),'$A_R$',str(replicateNumber),tempData[1]]) #change order here to change order of cells on tsplot
        preDf.append([float(i),'$B_R$',str(replicateNumber),tempData[0]])
        
        preDf.append([float(i),'Total',str(replicateNumber),tempData[2]])

odepreDf=[]
odepreDf2=[]

for k in range(cutoff):
    odepreDf.append([float(k),'ODE',replicates,data0[k]])
    
for k in range(cutoff):
    odepreDf2.append([float(k),'ODE2',replicates,data1[k]]) 
    
for k in range(cutoff):
    odepreDf2.append([float(k),'ODE3',replicates+1,totalODE[k]]) 

            
df=pandas.DataFrame(preDf,columns=['Time (days)', 'Cell Type','subject', 'Cell Number'])

odeDf=pandas.DataFrame(odepreDf,columns=['Time (days)', 'Cell Type','subject', 'Cell Number'])
odeDf2=pandas.DataFrame(odepreDf2,columns=['Time (days)', 'Cell Type','subject', 'Cell Number'])

#tsplot:

#sns.set_style("whitegrid", {'axes.grid' : False})

sns.set(font_scale=1.5)

ax = sns.tsplot(time="Time (days)", value="Cell Number",unit="subject", condition="Cell Type",data=df,err_style="unit_traces")

ax = sns.tsplot(time="Time (days)", value="Cell Number",unit="subject", condition="Cell Type",data=odeDf,color='k',linestyle='--',legend=True)
ax = sns.tsplot(time="Time (days)", value="Cell Number",unit="subject", condition="Cell Type",data=odeDf2,color='k',linestyle='--',legend=False)

#ax.set(yscale="log")

ax.grid(color='k',linestyle=':')

ax.set_xlim(0,200)
#ax.set_ylim(0,500)




#print(len(df))



#print(extinctionTimes)

data02=extinctionTimes

#######################################################################################################################################
#######################################################################################################################################
varFactor = 0.3

g_x_odd = -0.05 #Initial treatment (eg ALK TKI)
g_y_odd = 0.005

b_x_odd = (g_x_odd+varFactor)/2
b_y_odd = (g_y_odd+varFactor)/2

d_x_odd = (-g_x_odd+varFactor)/2
d_y_odd = (-g_y_odd+varFactor)/2

g_x_even = 0.005 
g_y_even = -0.05 

b_x_even = (g_x_even+varFactor)/2
b_y_even = (g_y_even+varFactor)/2

           
d_x_even = (-g_x_even+varFactor)/2
d_y_even = (-g_y_even+varFactor)/2

trans1 = 0.0001 
trans2 = trans1

relativePopInit = [9000,1000];

#Check for negative rates:
for parameter in [b_x_odd,b_y_odd,d_x_odd,d_y_odd,b_y_even,b_x_even,d_x_even,d_y_even]:
    if parameter < 0:
        raise NegativeRate

ApB=relativePopInit[1]/relativePopInit[0]
        
T_initial = np.log(((g_y_odd-g_x_odd)*(g_x_even-g_x_odd)+trans1*(g_y_odd+g_x_even-g_x_odd-g_y_even))/((trans1+ApB*(trans1+g_y_odd-g_x_odd))*(g_y_odd-g_y_even)))/(trans1+g_y_odd-g_x_odd)


k=(((g_y_odd-g_x_odd)*(g_x_even-g_x_odd)+trans1*(g_y_odd+g_x_even-g_x_odd-g_y_even))*(g_y_odd-g_y_even))/(((g_y_odd-g_y_even)*(g_x_even-g_y_even)+trans1*(g_y_odd+g_x_even-g_x_odd-g_y_even))*(g_x_even-g_x_odd))

T_treat_1 = 10.0
T_treat_2 = T_treat_1/k

T_treat = T_treat_1+T_treat_2

#T_initial -> T_treat_1 -> T_treat_2 -> T_treat_1

print(T_initial,T_treat_1,T_treat_2)



    





relativePop = relativePopInit


          
#-----------------Stoch----------------------------------------

def update(state):
    
    rt,re = random(), random()
    
    timeAfter = (state[2]-T_initial)/T_treat-int((state[2]-T_initial)/T_treat)
    
    
    if timeAfter>(T_treat_1/T_treat) or (state[2]<T_initial):
        a =(b_x_odd+d_x_odd+trans1)*state[0] + (b_y_odd+d_y_odd)*state[1]  
        
        if b_x_odd*state[0]/a>re:
            state[0] += 1
        elif (d_x_odd*state[0]+b_x_odd*state[0])/a>re:
            state[0] -= 1
        elif (b_y_odd*state[1]+d_x_odd*state[0]+b_x_odd*state[0])/a>re:
            state[1] += 1
        elif (d_y_odd*state[1]+b_y_odd*state[1]+d_x_odd*state[0]+b_x_odd*state[0])/a>re:
            state[1] -= 1
        else:
            state[0] -= 1
            state[1] += 1
    else:
        a =(b_x_even+d_x_even)*state[0] + (b_y_even+d_y_even+trans2)*state[1]
        
        if b_x_even*state[0]/a>re:
            state[0] += 1
        elif (d_x_even*state[0]+b_x_even*state[0])/a>re:
            state[0] -= 1
        elif (b_y_even*state[1]+d_x_even*state[0]+b_x_even*state[0])/a>re:
            state[1] += 1
        elif (d_y_even*state[1]+b_y_even*state[1]+d_x_even*state[0]+b_x_even*state[0])/a>re:
            state[1] -= 1
        else:
            state[0] += 1
            state[1] -= 1
             
    state[2] = state[2] - (1/a)*log(rt)
        
    
    return state



replicate = 0

extinctionTimes=[]

data=[]



while (replicate < replicates):   
    
    extinct=False;
    
    state = [relativePopInit[0], relativePopInit[1], 0.0] #initial conditions for state vector(0=sensitive, 1=resistant, 2=time)    
    
    life_history = []
    life_history.append(state)
    
    while (extinct==False) & (state[2]<cutoff):
    
        #print "Time: ", state[2]," Sensitive: ", state[0]," Resistant: ", state[1]
            
        new_state = copy(state)

        state = update(new_state)    
        
        life_history.append(state)
        
        if ((state[0]==0) & (state[1]==0)):
            extinct=True
            extinctionTimes.append(state[2])
            
    
    life_history.append([0,0,cutoff])
    
    data.append(life_history)
    
    replicate += 1
    
#-----------------ODE----------------------------------------


data0=[];
data1=[];

timeStep=0.01;

for k in range(0, int(cutoff/timeStep)):
    
    if (float(k)*timeStep).is_integer():
        data0.append(relativePop[0]);
        data1.append(relativePop[1]);
        
    r01=trans1;
    r10=trans2;    
    
    
    timeAfter = (k*timeStep-T_initial)/T_treat-int((k*timeStep-T_initial)/T_treat)
    
    if timeAfter>(T_treat_1/T_treat) or (k*timeStep<T_initial):
        g0=g_x_odd;
        g1=g_y_odd;
        r01=trans1;
        r10=0.0;
    else:
        g0=g_x_even;
        g1=g_y_even;
        r01=0.0;
        r10=trans2;
    
    
        
    relativePop[0]=relativePop[0]+(relativePop[0]*g0+relativePop[1]*r10-relativePop[0]*r01)*timeStep;
    relativePop[1]=relativePop[1]+(relativePop[1]*g1+relativePop[0]*r01-relativePop[1]*r10)*timeStep;

totalODE=[]

for i,na in enumerate(data0):
    totalODE.append(data0[i]+data1[i])    
    

                   
#------------------Plots--------------------------

preDf=[]


#Bin the data into replicates:

for replicateNumber,replicateData in enumerate(data):
    
    binData=[[[] for _ in range(cutoff+1)],[[] for _ in range(cutoff+1)]]
    
    for i,repData in enumerate(replicateData):
        
        binData[0][int(replicateData[i][2])].append(int(replicateData[i][0]))
        binData[1][int(replicateData[i][2])].append(int(replicateData[i][1]))
        
    for i,binList in enumerate(binData[0]):
        
        tempData=[0,0,0]        
        
        if len(binData[0][i])==0:
            tempData[0]=binData[0][i-1][len(binData[0][i-1])-1]
            binData[0][i].append(binData[0][i-1][len(binData[0][i-1])-1])
        else:
            tempData[0]=np.mean(binData[0][i])
        if len(binData[1][i])==0:    
            tempData[1]=binData[1][i-1][len(binData[1][i-1])-1]
            binData[1][i].append(binData[1][i-1][len(binData[0][i-1])-1])
        else:
            tempData[1]=np.mean(binData[1][i])
        
        tempData[2]=sum(tempData)
        
        preDf.append([float(i),'$A_R$',str(replicateNumber),tempData[1]]) #change order here to change order of cells on tsplot
        preDf.append([float(i),'$B_R$',str(replicateNumber),tempData[0]])
        
        preDf.append([float(i),'Total',str(replicateNumber),tempData[2]])

odepreDf=[]
odepreDf2=[]

for k in range(cutoff):
    odepreDf.append([float(k),'ODE',replicates,data0[k]])
    
for k in range(cutoff):
    odepreDf2.append([float(k),'ODE2',replicates,data1[k]]) 
    
for k in range(cutoff):
    odepreDf2.append([float(k),'ODE3',replicates+1,totalODE[k]]) 

            
df=pandas.DataFrame(preDf,columns=['Time (days)', 'Cell Type','subject', 'Cell Number'])

odeDf=pandas.DataFrame(odepreDf,columns=['Time (days)', 'Cell Type','subject', 'Cell Number'])
odeDf2=pandas.DataFrame(odepreDf2,columns=['Time (days)', 'Cell Type','subject', 'Cell Number'])

#tsplot:

#sns.set_style("whitegrid", {'axes.grid' : False})

sns.set(font_scale=1.5)

ax = sns.tsplot(time="Time (days)", value="Cell Number",unit="subject", condition="Cell Type",data=df,err_style="unit_traces")

ax = sns.tsplot(time="Time (days)", value="Cell Number",unit="subject", condition="Cell Type",data=odeDf,color='k',linestyle='--',legend=True)
ax = sns.tsplot(time="Time (days)", value="Cell Number",unit="subject", condition="Cell Type",data=odeDf2,color='k',linestyle='--',legend=False)

#ax.set(yscale="log")

ax.grid(color='k',linestyle=':')

ax.set_xlim(0,200)
#ax.set_ylim(0,500)




#print(len(df))



#print(extinctionTimes)

data03=extinctionTimes

#######################################################################################################################################
#######################################################################################################################################
varFactor = 0.4

g_x_odd = -0.05 #Initial treatment (eg ALK TKI)
g_y_odd = 0.005

b_x_odd = (g_x_odd+varFactor)/2
b_y_odd = (g_y_odd+varFactor)/2

d_x_odd = (-g_x_odd+varFactor)/2
d_y_odd = (-g_y_odd+varFactor)/2

g_x_even = 0.005 
g_y_even = -0.05 

b_x_even = (g_x_even+varFactor)/2
b_y_even = (g_y_even+varFactor)/2

           
d_x_even = (-g_x_even+varFactor)/2
d_y_even = (-g_y_even+varFactor)/2

trans1 = 0.0001 
trans2 = trans1

relativePopInit = [9000,1000];

#Check for negative rates:
for parameter in [b_x_odd,b_y_odd,d_x_odd,d_y_odd,b_y_even,b_x_even,d_x_even,d_y_even]:
    if parameter < 0:
        raise NegativeRate

ApB=relativePopInit[1]/relativePopInit[0]
        
T_initial = np.log(((g_y_odd-g_x_odd)*(g_x_even-g_x_odd)+trans1*(g_y_odd+g_x_even-g_x_odd-g_y_even))/((trans1+ApB*(trans1+g_y_odd-g_x_odd))*(g_y_odd-g_y_even)))/(trans1+g_y_odd-g_x_odd)


k=(((g_y_odd-g_x_odd)*(g_x_even-g_x_odd)+trans1*(g_y_odd+g_x_even-g_x_odd-g_y_even))*(g_y_odd-g_y_even))/(((g_y_odd-g_y_even)*(g_x_even-g_y_even)+trans1*(g_y_odd+g_x_even-g_x_odd-g_y_even))*(g_x_even-g_x_odd))

T_treat_1 = 10.0
T_treat_2 = T_treat_1/k

T_treat = T_treat_1+T_treat_2

#T_initial -> T_treat_1 -> T_treat_2 -> T_treat_1

print(T_initial,T_treat_1,T_treat_2)



    





relativePop = relativePopInit


          
#-----------------Stoch----------------------------------------

def update(state):
    
    rt,re = random(), random()
    
    timeAfter = (state[2]-T_initial)/T_treat-int((state[2]-T_initial)/T_treat)
    
    
    if timeAfter>(T_treat_1/T_treat) or (state[2]<T_initial):
        a =(b_x_odd+d_x_odd+trans1)*state[0] + (b_y_odd+d_y_odd)*state[1]  
        
        if b_x_odd*state[0]/a>re:
            state[0] += 1
        elif (d_x_odd*state[0]+b_x_odd*state[0])/a>re:
            state[0] -= 1
        elif (b_y_odd*state[1]+d_x_odd*state[0]+b_x_odd*state[0])/a>re:
            state[1] += 1
        elif (d_y_odd*state[1]+b_y_odd*state[1]+d_x_odd*state[0]+b_x_odd*state[0])/a>re:
            state[1] -= 1
        else:
            state[0] -= 1
            state[1] += 1
    else:
        a =(b_x_even+d_x_even)*state[0] + (b_y_even+d_y_even+trans2)*state[1]
        
        if b_x_even*state[0]/a>re:
            state[0] += 1
        elif (d_x_even*state[0]+b_x_even*state[0])/a>re:
            state[0] -= 1
        elif (b_y_even*state[1]+d_x_even*state[0]+b_x_even*state[0])/a>re:
            state[1] += 1
        elif (d_y_even*state[1]+b_y_even*state[1]+d_x_even*state[0]+b_x_even*state[0])/a>re:
            state[1] -= 1
        else:
            state[0] += 1
            state[1] -= 1
             
    state[2] = state[2] - (1/a)*log(rt)
        
    
    return state



replicate = 0

extinctionTimes=[]

data=[]



while (replicate < replicates):   
    
    extinct=False;
    
    state = [relativePopInit[0], relativePopInit[1], 0.0] #initial conditions for state vector(0=sensitive, 1=resistant, 2=time)    
    
    life_history = []
    life_history.append(state)
    
    while (extinct==False) & (state[2]<cutoff):
    
        #print "Time: ", state[2]," Sensitive: ", state[0]," Resistant: ", state[1]
            
        new_state = copy(state)

        state = update(new_state)    
        
        life_history.append(state)
        
        if ((state[0]==0) & (state[1]==0)):
            extinct=True
            extinctionTimes.append(state[2])
            
    
    life_history.append([0,0,cutoff])
    
    data.append(life_history)
    
    replicate += 1
    
#-----------------ODE----------------------------------------


data0=[];
data1=[];

timeStep=0.01;

for k in range(0, int(cutoff/timeStep)):
    
    if (float(k)*timeStep).is_integer():
        data0.append(relativePop[0]);
        data1.append(relativePop[1]);
        
    r01=trans1;
    r10=trans2;    
    
    
    timeAfter = (k*timeStep-T_initial)/T_treat-int((k*timeStep-T_initial)/T_treat)
    
    if timeAfter>(T_treat_1/T_treat) or (k*timeStep<T_initial):
        g0=g_x_odd;
        g1=g_y_odd;
        r01=trans1;
        r10=0.0;
    else:
        g0=g_x_even;
        g1=g_y_even;
        r01=0.0;
        r10=trans2;
    
    
        
    relativePop[0]=relativePop[0]+(relativePop[0]*g0+relativePop[1]*r10-relativePop[0]*r01)*timeStep;
    relativePop[1]=relativePop[1]+(relativePop[1]*g1+relativePop[0]*r01-relativePop[1]*r10)*timeStep;

totalODE=[]

for i,na in enumerate(data0):
    totalODE.append(data0[i]+data1[i])    
    

                   
#------------------Plots--------------------------

preDf=[]


#Bin the data into replicates:

for replicateNumber,replicateData in enumerate(data):
    
    binData=[[[] for _ in range(cutoff+1)],[[] for _ in range(cutoff+1)]]
    
    for i,repData in enumerate(replicateData):
        
        binData[0][int(replicateData[i][2])].append(int(replicateData[i][0]))
        binData[1][int(replicateData[i][2])].append(int(replicateData[i][1]))
        
    for i,binList in enumerate(binData[0]):
        
        tempData=[0,0,0]        
        
        if len(binData[0][i])==0:
            tempData[0]=binData[0][i-1][len(binData[0][i-1])-1]
            binData[0][i].append(binData[0][i-1][len(binData[0][i-1])-1])
        else:
            tempData[0]=np.mean(binData[0][i])
        if len(binData[1][i])==0:    
            tempData[1]=binData[1][i-1][len(binData[1][i-1])-1]
            binData[1][i].append(binData[1][i-1][len(binData[0][i-1])-1])
        else:
            tempData[1]=np.mean(binData[1][i])
        
        tempData[2]=sum(tempData)
        
        preDf.append([float(i),'$A_R$',str(replicateNumber),tempData[1]]) #change order here to change order of cells on tsplot
        preDf.append([float(i),'$B_R$',str(replicateNumber),tempData[0]])
        
        preDf.append([float(i),'Total',str(replicateNumber),tempData[2]])

odepreDf=[]
odepreDf2=[]

for k in range(cutoff):
    odepreDf.append([float(k),'ODE',replicates,data0[k]])
    
for k in range(cutoff):
    odepreDf2.append([float(k),'ODE2',replicates,data1[k]]) 
    
for k in range(cutoff):
    odepreDf2.append([float(k),'ODE3',replicates+1,totalODE[k]]) 

            
df=pandas.DataFrame(preDf,columns=['Time (days)', 'Cell Type','subject', 'Cell Number'])

odeDf=pandas.DataFrame(odepreDf,columns=['Time (days)', 'Cell Type','subject', 'Cell Number'])
odeDf2=pandas.DataFrame(odepreDf2,columns=['Time (days)', 'Cell Type','subject', 'Cell Number'])

#tsplot:

#sns.set_style("whitegrid", {'axes.grid' : False})

sns.set(font_scale=1.5)

ax = sns.tsplot(time="Time (days)", value="Cell Number",unit="subject", condition="Cell Type",data=df,err_style="unit_traces")

ax = sns.tsplot(time="Time (days)", value="Cell Number",unit="subject", condition="Cell Type",data=odeDf,color='k',linestyle='--',legend=True)
ax = sns.tsplot(time="Time (days)", value="Cell Number",unit="subject", condition="Cell Type",data=odeDf2,color='k',linestyle='--',legend=False)

#ax.set(yscale="log")

ax.grid(color='k',linestyle=':')

ax.set_xlim(0,200)
#ax.set_ylim(0,500)




#print(len(df))



#print(extinctionTimes)

data04=extinctionTimes

#######################################################################################################################################
#######################################################################################################################################
varFactor = 0.5

g_x_odd = -0.05 #Initial treatment (eg ALK TKI)
g_y_odd = 0.005

b_x_odd = (g_x_odd+varFactor)/2
b_y_odd = (g_y_odd+varFactor)/2

d_x_odd = (-g_x_odd+varFactor)/2
d_y_odd = (-g_y_odd+varFactor)/2

g_x_even = 0.005 
g_y_even = -0.05 

b_x_even = (g_x_even+varFactor)/2
b_y_even = (g_y_even+varFactor)/2

           
d_x_even = (-g_x_even+varFactor)/2
d_y_even = (-g_y_even+varFactor)/2

trans1 = 0.0001 
trans2 = trans1

relativePopInit = [9000,1000];

#Check for negative rates:
for parameter in [b_x_odd,b_y_odd,d_x_odd,d_y_odd,b_y_even,b_x_even,d_x_even,d_y_even]:
    if parameter < 0:
        raise NegativeRate

ApB=relativePopInit[1]/relativePopInit[0]
        
T_initial = np.log(((g_y_odd-g_x_odd)*(g_x_even-g_x_odd)+trans1*(g_y_odd+g_x_even-g_x_odd-g_y_even))/((trans1+ApB*(trans1+g_y_odd-g_x_odd))*(g_y_odd-g_y_even)))/(trans1+g_y_odd-g_x_odd)


k=(((g_y_odd-g_x_odd)*(g_x_even-g_x_odd)+trans1*(g_y_odd+g_x_even-g_x_odd-g_y_even))*(g_y_odd-g_y_even))/(((g_y_odd-g_y_even)*(g_x_even-g_y_even)+trans1*(g_y_odd+g_x_even-g_x_odd-g_y_even))*(g_x_even-g_x_odd))

T_treat_1 = 10.0
T_treat_2 = T_treat_1/k

T_treat = T_treat_1+T_treat_2

#T_initial -> T_treat_1 -> T_treat_2 -> T_treat_1

print(T_initial,T_treat_1,T_treat_2)



    





relativePop = relativePopInit


          
#-----------------Stoch----------------------------------------

def update(state):
    
    rt,re = random(), random()
    
    timeAfter = (state[2]-T_initial)/T_treat-int((state[2]-T_initial)/T_treat)
    
    
    if timeAfter>(T_treat_1/T_treat) or (state[2]<T_initial):
        a =(b_x_odd+d_x_odd+trans1)*state[0] + (b_y_odd+d_y_odd)*state[1]  
        
        if b_x_odd*state[0]/a>re:
            state[0] += 1
        elif (d_x_odd*state[0]+b_x_odd*state[0])/a>re:
            state[0] -= 1
        elif (b_y_odd*state[1]+d_x_odd*state[0]+b_x_odd*state[0])/a>re:
            state[1] += 1
        elif (d_y_odd*state[1]+b_y_odd*state[1]+d_x_odd*state[0]+b_x_odd*state[0])/a>re:
            state[1] -= 1
        else:
            state[0] -= 1
            state[1] += 1
    else:
        a =(b_x_even+d_x_even)*state[0] + (b_y_even+d_y_even+trans2)*state[1]
        
        if b_x_even*state[0]/a>re:
            state[0] += 1
        elif (d_x_even*state[0]+b_x_even*state[0])/a>re:
            state[0] -= 1
        elif (b_y_even*state[1]+d_x_even*state[0]+b_x_even*state[0])/a>re:
            state[1] += 1
        elif (d_y_even*state[1]+b_y_even*state[1]+d_x_even*state[0]+b_x_even*state[0])/a>re:
            state[1] -= 1
        else:
            state[0] += 1
            state[1] -= 1
             
    state[2] = state[2] - (1/a)*log(rt)
        
    
    return state



replicate = 0

extinctionTimes=[]

data=[]



while (replicate < replicates):   
    
    extinct=False;
    
    state = [relativePopInit[0], relativePopInit[1], 0.0] #initial conditions for state vector(0=sensitive, 1=resistant, 2=time)    
    
    life_history = []
    life_history.append(state)
    
    while (extinct==False) & (state[2]<cutoff):
    
        #print "Time: ", state[2]," Sensitive: ", state[0]," Resistant: ", state[1]
            
        new_state = copy(state)

        state = update(new_state)    
        
        life_history.append(state)
        
        if ((state[0]==0) & (state[1]==0)):
            extinct=True
            extinctionTimes.append(state[2])
            
    
    life_history.append([0,0,cutoff])
    
    data.append(life_history)
    
    replicate += 1
    
#-----------------ODE----------------------------------------


data0=[];
data1=[];

timeStep=0.01;

for k in range(0, int(cutoff/timeStep)):
    
    if (float(k)*timeStep).is_integer():
        data0.append(relativePop[0]);
        data1.append(relativePop[1]);
        
    r01=trans1;
    r10=trans2;    
    
    
    timeAfter = (k*timeStep-T_initial)/T_treat-int((k*timeStep-T_initial)/T_treat)
    
    if timeAfter>(T_treat_1/T_treat) or (k*timeStep<T_initial):
        g0=g_x_odd;
        g1=g_y_odd;
        r01=trans1;
        r10=0.0;
    else:
        g0=g_x_even;
        g1=g_y_even;
        r01=0.0;
        r10=trans2;
    
    
        
    relativePop[0]=relativePop[0]+(relativePop[0]*g0+relativePop[1]*r10-relativePop[0]*r01)*timeStep;
    relativePop[1]=relativePop[1]+(relativePop[1]*g1+relativePop[0]*r01-relativePop[1]*r10)*timeStep;

totalODE=[]

for i,na in enumerate(data0):
    totalODE.append(data0[i]+data1[i])    
    

                   
#------------------Plots--------------------------

preDf=[]


#Bin the data into replicates:

for replicateNumber,replicateData in enumerate(data):
    
    binData=[[[] for _ in range(cutoff+1)],[[] for _ in range(cutoff+1)]]
    
    for i,repData in enumerate(replicateData):
        
        binData[0][int(replicateData[i][2])].append(int(replicateData[i][0]))
        binData[1][int(replicateData[i][2])].append(int(replicateData[i][1]))
        
    for i,binList in enumerate(binData[0]):
        
        tempData=[0,0,0]        
        
        if len(binData[0][i])==0:
            tempData[0]=binData[0][i-1][len(binData[0][i-1])-1]
            binData[0][i].append(binData[0][i-1][len(binData[0][i-1])-1])
        else:
            tempData[0]=np.mean(binData[0][i])
        if len(binData[1][i])==0:    
            tempData[1]=binData[1][i-1][len(binData[1][i-1])-1]
            binData[1][i].append(binData[1][i-1][len(binData[0][i-1])-1])
        else:
            tempData[1]=np.mean(binData[1][i])
        
        tempData[2]=sum(tempData)
        
        preDf.append([float(i),'$A_R$',str(replicateNumber),tempData[1]]) #change order here to change order of cells on tsplot
        preDf.append([float(i),'$B_R$',str(replicateNumber),tempData[0]])
        
        preDf.append([float(i),'Total',str(replicateNumber),tempData[2]])

odepreDf=[]
odepreDf2=[]

for k in range(cutoff):
    odepreDf.append([float(k),'ODE',replicates,data0[k]])
    
for k in range(cutoff):
    odepreDf2.append([float(k),'ODE2',replicates,data1[k]]) 
    
for k in range(cutoff):
    odepreDf2.append([float(k),'ODE3',replicates+1,totalODE[k]]) 

            
df=pandas.DataFrame(preDf,columns=['Time (days)', 'Cell Type','subject', 'Cell Number'])

odeDf=pandas.DataFrame(odepreDf,columns=['Time (days)', 'Cell Type','subject', 'Cell Number'])
odeDf2=pandas.DataFrame(odepreDf2,columns=['Time (days)', 'Cell Type','subject', 'Cell Number'])

#tsplot:

#sns.set_style("whitegrid", {'axes.grid' : False})

sns.set(font_scale=1.5)

ax = sns.tsplot(time="Time (days)", value="Cell Number",unit="subject", condition="Cell Type",data=df,err_style="unit_traces")

ax = sns.tsplot(time="Time (days)", value="Cell Number",unit="subject", condition="Cell Type",data=odeDf,color='k',linestyle='--',legend=True)
ax = sns.tsplot(time="Time (days)", value="Cell Number",unit="subject", condition="Cell Type",data=odeDf2,color='k',linestyle='--',legend=False)

#ax.set(yscale="log")

ax.grid(color='k',linestyle=':')

ax.set_xlim(0,200)
#ax.set_ylim(0,500)




#print(len(df))



#print(extinctionTimes)

data05=extinctionTimes

#######################################################################################################################################
#######################################################################################################################################
varFactor = 0.6

g_x_odd = -0.05 #Initial treatment (eg ALK TKI)
g_y_odd = 0.005

b_x_odd = (g_x_odd+varFactor)/2
b_y_odd = (g_y_odd+varFactor)/2

d_x_odd = (-g_x_odd+varFactor)/2
d_y_odd = (-g_y_odd+varFactor)/2

g_x_even = 0.005 
g_y_even = -0.05 

b_x_even = (g_x_even+varFactor)/2
b_y_even = (g_y_even+varFactor)/2

           
d_x_even = (-g_x_even+varFactor)/2
d_y_even = (-g_y_even+varFactor)/2

trans1 = 0.0001 
trans2 = trans1

relativePopInit = [9000,1000];

#Check for negative rates:
for parameter in [b_x_odd,b_y_odd,d_x_odd,d_y_odd,b_y_even,b_x_even,d_x_even,d_y_even]:
    if parameter < 0:
        raise NegativeRate

ApB=relativePopInit[1]/relativePopInit[0]
        
T_initial = np.log(((g_y_odd-g_x_odd)*(g_x_even-g_x_odd)+trans1*(g_y_odd+g_x_even-g_x_odd-g_y_even))/((trans1+ApB*(trans1+g_y_odd-g_x_odd))*(g_y_odd-g_y_even)))/(trans1+g_y_odd-g_x_odd)


k=(((g_y_odd-g_x_odd)*(g_x_even-g_x_odd)+trans1*(g_y_odd+g_x_even-g_x_odd-g_y_even))*(g_y_odd-g_y_even))/(((g_y_odd-g_y_even)*(g_x_even-g_y_even)+trans1*(g_y_odd+g_x_even-g_x_odd-g_y_even))*(g_x_even-g_x_odd))

T_treat_1 = 10.0
T_treat_2 = T_treat_1/k

T_treat = T_treat_1+T_treat_2

#T_initial -> T_treat_1 -> T_treat_2 -> T_treat_1

print(T_initial,T_treat_1,T_treat_2)



    





relativePop = relativePopInit


          
#-----------------Stoch----------------------------------------

def update(state):
    
    rt,re = random(), random()
    
    timeAfter = (state[2]-T_initial)/T_treat-int((state[2]-T_initial)/T_treat)
    
    
    if timeAfter>(T_treat_1/T_treat) or (state[2]<T_initial):
        a =(b_x_odd+d_x_odd+trans1)*state[0] + (b_y_odd+d_y_odd)*state[1]  
        
        if b_x_odd*state[0]/a>re:
            state[0] += 1
        elif (d_x_odd*state[0]+b_x_odd*state[0])/a>re:
            state[0] -= 1
        elif (b_y_odd*state[1]+d_x_odd*state[0]+b_x_odd*state[0])/a>re:
            state[1] += 1
        elif (d_y_odd*state[1]+b_y_odd*state[1]+d_x_odd*state[0]+b_x_odd*state[0])/a>re:
            state[1] -= 1
        else:
            state[0] -= 1
            state[1] += 1
    else:
        a =(b_x_even+d_x_even)*state[0] + (b_y_even+d_y_even+trans2)*state[1]
        
        if b_x_even*state[0]/a>re:
            state[0] += 1
        elif (d_x_even*state[0]+b_x_even*state[0])/a>re:
            state[0] -= 1
        elif (b_y_even*state[1]+d_x_even*state[0]+b_x_even*state[0])/a>re:
            state[1] += 1
        elif (d_y_even*state[1]+b_y_even*state[1]+d_x_even*state[0]+b_x_even*state[0])/a>re:
            state[1] -= 1
        else:
            state[0] += 1
            state[1] -= 1
             
    state[2] = state[2] - (1/a)*log(rt)
        
    
    return state



replicate = 0

extinctionTimes=[]

data=[]



while (replicate < replicates):   
    
    extinct=False;
    
    state = [relativePopInit[0], relativePopInit[1], 0.0] #initial conditions for state vector(0=sensitive, 1=resistant, 2=time)    
    
    life_history = []
    life_history.append(state)
    
    while (extinct==False) & (state[2]<cutoff):
    
        #print "Time: ", state[2]," Sensitive: ", state[0]," Resistant: ", state[1]
            
        new_state = copy(state)

        state = update(new_state)    
        
        life_history.append(state)
        
        if ((state[0]==0) & (state[1]==0)):
            extinct=True
            extinctionTimes.append(state[2])
            
    
    life_history.append([0,0,cutoff])
    
    data.append(life_history)
    
    replicate += 1
    
#-----------------ODE----------------------------------------


data0=[];
data1=[];

timeStep=0.01;

for k in range(0, int(cutoff/timeStep)):
    
    if (float(k)*timeStep).is_integer():
        data0.append(relativePop[0]);
        data1.append(relativePop[1]);
        
    r01=trans1;
    r10=trans2;    
    
    
    timeAfter = (k*timeStep-T_initial)/T_treat-int((k*timeStep-T_initial)/T_treat)
    
    if timeAfter>(T_treat_1/T_treat) or (k*timeStep<T_initial):
        g0=g_x_odd;
        g1=g_y_odd;
        r01=trans1;
        r10=0.0;
    else:
        g0=g_x_even;
        g1=g_y_even;
        r01=0.0;
        r10=trans2;
    
    
        
    relativePop[0]=relativePop[0]+(relativePop[0]*g0+relativePop[1]*r10-relativePop[0]*r01)*timeStep;
    relativePop[1]=relativePop[1]+(relativePop[1]*g1+relativePop[0]*r01-relativePop[1]*r10)*timeStep;

totalODE=[]

for i,na in enumerate(data0):
    totalODE.append(data0[i]+data1[i])    
    

                   
#------------------Plots--------------------------

preDf=[]


#Bin the data into replicates:

for replicateNumber,replicateData in enumerate(data):
    
    binData=[[[] for _ in range(cutoff+1)],[[] for _ in range(cutoff+1)]]
    
    for i,repData in enumerate(replicateData):
        
        binData[0][int(replicateData[i][2])].append(int(replicateData[i][0]))
        binData[1][int(replicateData[i][2])].append(int(replicateData[i][1]))
        
    for i,binList in enumerate(binData[0]):
        
        tempData=[0,0,0]        
        
        if len(binData[0][i])==0:
            tempData[0]=binData[0][i-1][len(binData[0][i-1])-1]
            binData[0][i].append(binData[0][i-1][len(binData[0][i-1])-1])
        else:
            tempData[0]=np.mean(binData[0][i])
        if len(binData[1][i])==0:    
            tempData[1]=binData[1][i-1][len(binData[1][i-1])-1]
            binData[1][i].append(binData[1][i-1][len(binData[0][i-1])-1])
        else:
            tempData[1]=np.mean(binData[1][i])
        
        tempData[2]=sum(tempData)
        
        preDf.append([float(i),'$A_R$',str(replicateNumber),tempData[1]]) #change order here to change order of cells on tsplot
        preDf.append([float(i),'$B_R$',str(replicateNumber),tempData[0]])
        
        preDf.append([float(i),'Total',str(replicateNumber),tempData[2]])

odepreDf=[]
odepreDf2=[]

for k in range(cutoff):
    odepreDf.append([float(k),'ODE',replicates,data0[k]])
    
for k in range(cutoff):
    odepreDf2.append([float(k),'ODE2',replicates,data1[k]]) 
    
for k in range(cutoff):
    odepreDf2.append([float(k),'ODE3',replicates+1,totalODE[k]]) 

            
df=pandas.DataFrame(preDf,columns=['Time (days)', 'Cell Type','subject', 'Cell Number'])

odeDf=pandas.DataFrame(odepreDf,columns=['Time (days)', 'Cell Type','subject', 'Cell Number'])
odeDf2=pandas.DataFrame(odepreDf2,columns=['Time (days)', 'Cell Type','subject', 'Cell Number'])

#tsplot:

#sns.set_style("whitegrid", {'axes.grid' : False})

sns.set(font_scale=1.5)

ax = sns.tsplot(time="Time (days)", value="Cell Number",unit="subject", condition="Cell Type",data=df,err_style="unit_traces")

ax = sns.tsplot(time="Time (days)", value="Cell Number",unit="subject", condition="Cell Type",data=odeDf,color='k',linestyle='--',legend=True)
ax = sns.tsplot(time="Time (days)", value="Cell Number",unit="subject", condition="Cell Type",data=odeDf2,color='k',linestyle='--',legend=False)

#ax.set(yscale="log")

ax.grid(color='k',linestyle=':')

ax.set_xlim(0,200)
#ax.set_ylim(0,500)




#print(len(df))



#print(extinctionTimes)

data06=extinctionTimes

#######################################################################################################################################
#######################################################################################################################################
varFactor = 0.7

g_x_odd = -0.05 #Initial treatment (eg ALK TKI)
g_y_odd = 0.005

b_x_odd = (g_x_odd+varFactor)/2
b_y_odd = (g_y_odd+varFactor)/2

d_x_odd = (-g_x_odd+varFactor)/2
d_y_odd = (-g_y_odd+varFactor)/2

g_x_even = 0.005 
g_y_even = -0.05 

b_x_even = (g_x_even+varFactor)/2
b_y_even = (g_y_even+varFactor)/2

           
d_x_even = (-g_x_even+varFactor)/2
d_y_even = (-g_y_even+varFactor)/2

trans1 = 0.0001 
trans2 = trans1

relativePopInit = [9000,1000];

#Check for negative rates:
for parameter in [b_x_odd,b_y_odd,d_x_odd,d_y_odd,b_y_even,b_x_even,d_x_even,d_y_even]:
    if parameter < 0:
        raise NegativeRate

ApB=relativePopInit[1]/relativePopInit[0]
        
T_initial = np.log(((g_y_odd-g_x_odd)*(g_x_even-g_x_odd)+trans1*(g_y_odd+g_x_even-g_x_odd-g_y_even))/((trans1+ApB*(trans1+g_y_odd-g_x_odd))*(g_y_odd-g_y_even)))/(trans1+g_y_odd-g_x_odd)


k=(((g_y_odd-g_x_odd)*(g_x_even-g_x_odd)+trans1*(g_y_odd+g_x_even-g_x_odd-g_y_even))*(g_y_odd-g_y_even))/(((g_y_odd-g_y_even)*(g_x_even-g_y_even)+trans1*(g_y_odd+g_x_even-g_x_odd-g_y_even))*(g_x_even-g_x_odd))

T_treat_1 = 10.0
T_treat_2 = T_treat_1/k

T_treat = T_treat_1+T_treat_2

#T_initial -> T_treat_1 -> T_treat_2 -> T_treat_1

print(T_initial,T_treat_1,T_treat_2)



    





relativePop = relativePopInit


          
#-----------------Stoch----------------------------------------

def update(state):
    
    rt,re = random(), random()
    
    timeAfter = (state[2]-T_initial)/T_treat-int((state[2]-T_initial)/T_treat)
    
    
    if timeAfter>(T_treat_1/T_treat) or (state[2]<T_initial):
        a =(b_x_odd+d_x_odd+trans1)*state[0] + (b_y_odd+d_y_odd)*state[1]  
        
        if b_x_odd*state[0]/a>re:
            state[0] += 1
        elif (d_x_odd*state[0]+b_x_odd*state[0])/a>re:
            state[0] -= 1
        elif (b_y_odd*state[1]+d_x_odd*state[0]+b_x_odd*state[0])/a>re:
            state[1] += 1
        elif (d_y_odd*state[1]+b_y_odd*state[1]+d_x_odd*state[0]+b_x_odd*state[0])/a>re:
            state[1] -= 1
        else:
            state[0] -= 1
            state[1] += 1
    else:
        a =(b_x_even+d_x_even)*state[0] + (b_y_even+d_y_even+trans2)*state[1]
        
        if b_x_even*state[0]/a>re:
            state[0] += 1
        elif (d_x_even*state[0]+b_x_even*state[0])/a>re:
            state[0] -= 1
        elif (b_y_even*state[1]+d_x_even*state[0]+b_x_even*state[0])/a>re:
            state[1] += 1
        elif (d_y_even*state[1]+b_y_even*state[1]+d_x_even*state[0]+b_x_even*state[0])/a>re:
            state[1] -= 1
        else:
            state[0] += 1
            state[1] -= 1
             
    state[2] = state[2] - (1/a)*log(rt)
        
    
    return state



replicate = 0

extinctionTimes=[]

data=[]



while (replicate < replicates):   
    
    extinct=False;
    
    state = [relativePopInit[0], relativePopInit[1], 0.0] #initial conditions for state vector(0=sensitive, 1=resistant, 2=time)    
    
    life_history = []
    life_history.append(state)
    
    while (extinct==False) & (state[2]<cutoff):
    
        #print "Time: ", state[2]," Sensitive: ", state[0]," Resistant: ", state[1]
            
        new_state = copy(state)

        state = update(new_state)    
        
        life_history.append(state)
        
        if ((state[0]==0) & (state[1]==0)):
            extinct=True
            extinctionTimes.append(state[2])
            
    
    life_history.append([0,0,cutoff])
    
    data.append(life_history)
    
    replicate += 1
    
#-----------------ODE----------------------------------------


data0=[];
data1=[];

timeStep=0.01;

for k in range(0, int(cutoff/timeStep)):
    
    if (float(k)*timeStep).is_integer():
        data0.append(relativePop[0]);
        data1.append(relativePop[1]);
        
    r01=trans1;
    r10=trans2;    
    
    
    timeAfter = (k*timeStep-T_initial)/T_treat-int((k*timeStep-T_initial)/T_treat)
    
    if timeAfter>(T_treat_1/T_treat) or (k*timeStep<T_initial):
        g0=g_x_odd;
        g1=g_y_odd;
        r01=trans1;
        r10=0.0;
    else:
        g0=g_x_even;
        g1=g_y_even;
        r01=0.0;
        r10=trans2;
    
    
        
    relativePop[0]=relativePop[0]+(relativePop[0]*g0+relativePop[1]*r10-relativePop[0]*r01)*timeStep;
    relativePop[1]=relativePop[1]+(relativePop[1]*g1+relativePop[0]*r01-relativePop[1]*r10)*timeStep;

totalODE=[]

for i,na in enumerate(data0):
    totalODE.append(data0[i]+data1[i])    
    

                   
#------------------Plots--------------------------

preDf=[]


#Bin the data into replicates:

for replicateNumber,replicateData in enumerate(data):
    
    binData=[[[] for _ in range(cutoff+1)],[[] for _ in range(cutoff+1)]]
    
    for i,repData in enumerate(replicateData):
        
        binData[0][int(replicateData[i][2])].append(int(replicateData[i][0]))
        binData[1][int(replicateData[i][2])].append(int(replicateData[i][1]))
        
    for i,binList in enumerate(binData[0]):
        
        tempData=[0,0,0]        
        
        if len(binData[0][i])==0:
            tempData[0]=binData[0][i-1][len(binData[0][i-1])-1]
            binData[0][i].append(binData[0][i-1][len(binData[0][i-1])-1])
        else:
            tempData[0]=np.mean(binData[0][i])
        if len(binData[1][i])==0:    
            tempData[1]=binData[1][i-1][len(binData[1][i-1])-1]
            binData[1][i].append(binData[1][i-1][len(binData[0][i-1])-1])
        else:
            tempData[1]=np.mean(binData[1][i])
        
        tempData[2]=sum(tempData)
        
        preDf.append([float(i),'$A_R$',str(replicateNumber),tempData[1]]) #change order here to change order of cells on tsplot
        preDf.append([float(i),'$B_R$',str(replicateNumber),tempData[0]])
        
        preDf.append([float(i),'Total',str(replicateNumber),tempData[2]])

odepreDf=[]
odepreDf2=[]

for k in range(cutoff):
    odepreDf.append([float(k),'ODE',replicates,data0[k]])
    
for k in range(cutoff):
    odepreDf2.append([float(k),'ODE2',replicates,data1[k]]) 
    
for k in range(cutoff):
    odepreDf2.append([float(k),'ODE3',replicates+1,totalODE[k]]) 

            
df=pandas.DataFrame(preDf,columns=['Time (days)', 'Cell Type','subject', 'Cell Number'])

odeDf=pandas.DataFrame(odepreDf,columns=['Time (days)', 'Cell Type','subject', 'Cell Number'])
odeDf2=pandas.DataFrame(odepreDf2,columns=['Time (days)', 'Cell Type','subject', 'Cell Number'])

#tsplot:

#sns.set_style("whitegrid", {'axes.grid' : False})

sns.set(font_scale=1.5)

ax = sns.tsplot(time="Time (days)", value="Cell Number",unit="subject", condition="Cell Type",data=df,err_style="unit_traces")

ax = sns.tsplot(time="Time (days)", value="Cell Number",unit="subject", condition="Cell Type",data=odeDf,color='k',linestyle='--',legend=True)
ax = sns.tsplot(time="Time (days)", value="Cell Number",unit="subject", condition="Cell Type",data=odeDf2,color='k',linestyle='--',legend=False)

#ax.set(yscale="log")

ax.grid(color='k',linestyle=':')

ax.set_xlim(0,200)
#ax.set_ylim(0,500)




#print(len(df))



#print(extinctionTimes)

data07=extinctionTimes

#######################################################################################################################################
#######################################################################################################################################
varFactor = 0.8

g_x_odd = -0.05 #Initial treatment (eg ALK TKI)
g_y_odd = 0.005

b_x_odd = (g_x_odd+varFactor)/2
b_y_odd = (g_y_odd+varFactor)/2

d_x_odd = (-g_x_odd+varFactor)/2
d_y_odd = (-g_y_odd+varFactor)/2

g_x_even = 0.005 
g_y_even = -0.05 

b_x_even = (g_x_even+varFactor)/2
b_y_even = (g_y_even+varFactor)/2

           
d_x_even = (-g_x_even+varFactor)/2
d_y_even = (-g_y_even+varFactor)/2

trans1 = 0.0001 
trans2 = trans1

relativePopInit = [9000,1000];

#Check for negative rates:
for parameter in [b_x_odd,b_y_odd,d_x_odd,d_y_odd,b_y_even,b_x_even,d_x_even,d_y_even]:
    if parameter < 0:
        raise NegativeRate

ApB=relativePopInit[1]/relativePopInit[0]
        
T_initial = np.log(((g_y_odd-g_x_odd)*(g_x_even-g_x_odd)+trans1*(g_y_odd+g_x_even-g_x_odd-g_y_even))/((trans1+ApB*(trans1+g_y_odd-g_x_odd))*(g_y_odd-g_y_even)))/(trans1+g_y_odd-g_x_odd)


k=(((g_y_odd-g_x_odd)*(g_x_even-g_x_odd)+trans1*(g_y_odd+g_x_even-g_x_odd-g_y_even))*(g_y_odd-g_y_even))/(((g_y_odd-g_y_even)*(g_x_even-g_y_even)+trans1*(g_y_odd+g_x_even-g_x_odd-g_y_even))*(g_x_even-g_x_odd))

T_treat_1 = 10.0
T_treat_2 = T_treat_1/k

T_treat = T_treat_1+T_treat_2

#T_initial -> T_treat_1 -> T_treat_2 -> T_treat_1

print(T_initial,T_treat_1,T_treat_2)



    





relativePop = relativePopInit


          
#-----------------Stoch----------------------------------------

def update(state):
    
    rt,re = random(), random()
    
    timeAfter = (state[2]-T_initial)/T_treat-int((state[2]-T_initial)/T_treat)
    
    
    if timeAfter>(T_treat_1/T_treat) or (state[2]<T_initial):
        a =(b_x_odd+d_x_odd+trans1)*state[0] + (b_y_odd+d_y_odd)*state[1]  
        
        if b_x_odd*state[0]/a>re:
            state[0] += 1
        elif (d_x_odd*state[0]+b_x_odd*state[0])/a>re:
            state[0] -= 1
        elif (b_y_odd*state[1]+d_x_odd*state[0]+b_x_odd*state[0])/a>re:
            state[1] += 1
        elif (d_y_odd*state[1]+b_y_odd*state[1]+d_x_odd*state[0]+b_x_odd*state[0])/a>re:
            state[1] -= 1
        else:
            state[0] -= 1
            state[1] += 1
    else:
        a =(b_x_even+d_x_even)*state[0] + (b_y_even+d_y_even+trans2)*state[1]
        
        if b_x_even*state[0]/a>re:
            state[0] += 1
        elif (d_x_even*state[0]+b_x_even*state[0])/a>re:
            state[0] -= 1
        elif (b_y_even*state[1]+d_x_even*state[0]+b_x_even*state[0])/a>re:
            state[1] += 1
        elif (d_y_even*state[1]+b_y_even*state[1]+d_x_even*state[0]+b_x_even*state[0])/a>re:
            state[1] -= 1
        else:
            state[0] += 1
            state[1] -= 1
             
    state[2] = state[2] - (1/a)*log(rt)
        
    
    return state



replicate = 0

extinctionTimes=[]

data=[]



while (replicate < replicates):   
    
    extinct=False;
    
    state = [relativePopInit[0], relativePopInit[1], 0.0] #initial conditions for state vector(0=sensitive, 1=resistant, 2=time)    
    
    life_history = []
    life_history.append(state)
    
    while (extinct==False) & (state[2]<cutoff):
    
        #print "Time: ", state[2]," Sensitive: ", state[0]," Resistant: ", state[1]
            
        new_state = copy(state)

        state = update(new_state)    
        
        life_history.append(state)
        
        if ((state[0]==0) & (state[1]==0)):
            extinct=True
            extinctionTimes.append(state[2])
            
    
    life_history.append([0,0,cutoff])
    
    data.append(life_history)
    
    replicate += 1
    
#-----------------ODE----------------------------------------


data0=[];
data1=[];

timeStep=0.01;

for k in range(0, int(cutoff/timeStep)):
    
    if (float(k)*timeStep).is_integer():
        data0.append(relativePop[0]);
        data1.append(relativePop[1]);
        
    r01=trans1;
    r10=trans2;    
    
    
    timeAfter = (k*timeStep-T_initial)/T_treat-int((k*timeStep-T_initial)/T_treat)
    
    if timeAfter>(T_treat_1/T_treat) or (k*timeStep<T_initial):
        g0=g_x_odd;
        g1=g_y_odd;
        r01=trans1;
        r10=0.0;
    else:
        g0=g_x_even;
        g1=g_y_even;
        r01=0.0;
        r10=trans2;
    
    
        
    relativePop[0]=relativePop[0]+(relativePop[0]*g0+relativePop[1]*r10-relativePop[0]*r01)*timeStep;
    relativePop[1]=relativePop[1]+(relativePop[1]*g1+relativePop[0]*r01-relativePop[1]*r10)*timeStep;

totalODE=[]

for i,na in enumerate(data0):
    totalODE.append(data0[i]+data1[i])    
    

                   
#------------------Plots--------------------------

preDf=[]


#Bin the data into replicates:

for replicateNumber,replicateData in enumerate(data):
    
    binData=[[[] for _ in range(cutoff+1)],[[] for _ in range(cutoff+1)]]
    
    for i,repData in enumerate(replicateData):
        
        binData[0][int(replicateData[i][2])].append(int(replicateData[i][0]))
        binData[1][int(replicateData[i][2])].append(int(replicateData[i][1]))
        
    for i,binList in enumerate(binData[0]):
        
        tempData=[0,0,0]        
        
        if len(binData[0][i])==0:
            tempData[0]=binData[0][i-1][len(binData[0][i-1])-1]
            binData[0][i].append(binData[0][i-1][len(binData[0][i-1])-1])
        else:
            tempData[0]=np.mean(binData[0][i])
        if len(binData[1][i])==0:    
            tempData[1]=binData[1][i-1][len(binData[1][i-1])-1]
            binData[1][i].append(binData[1][i-1][len(binData[0][i-1])-1])
        else:
            tempData[1]=np.mean(binData[1][i])
        
        tempData[2]=sum(tempData)
        
        preDf.append([float(i),'$A_R$',str(replicateNumber),tempData[1]]) #change order here to change order of cells on tsplot
        preDf.append([float(i),'$B_R$',str(replicateNumber),tempData[0]])
        
        preDf.append([float(i),'Total',str(replicateNumber),tempData[2]])

odepreDf=[]
odepreDf2=[]

for k in range(cutoff):
    odepreDf.append([float(k),'ODE',replicates,data0[k]])
    
for k in range(cutoff):
    odepreDf2.append([float(k),'ODE2',replicates,data1[k]]) 
    
for k in range(cutoff):
    odepreDf2.append([float(k),'ODE3',replicates+1,totalODE[k]]) 

            
df=pandas.DataFrame(preDf,columns=['Time (days)', 'Cell Type','subject', 'Cell Number'])

odeDf=pandas.DataFrame(odepreDf,columns=['Time (days)', 'Cell Type','subject', 'Cell Number'])
odeDf2=pandas.DataFrame(odepreDf2,columns=['Time (days)', 'Cell Type','subject', 'Cell Number'])

#tsplot:

#sns.set_style("whitegrid", {'axes.grid' : False})

sns.set(font_scale=1.5)

ax = sns.tsplot(time="Time (days)", value="Cell Number",unit="subject", condition="Cell Type",data=df,err_style="unit_traces")

ax = sns.tsplot(time="Time (days)", value="Cell Number",unit="subject", condition="Cell Type",data=odeDf,color='k',linestyle='--',legend=True)
ax = sns.tsplot(time="Time (days)", value="Cell Number",unit="subject", condition="Cell Type",data=odeDf2,color='k',linestyle='--',legend=False)

#ax.set(yscale="log")

ax.grid(color='k',linestyle=':')

ax.set_xlim(0,200)
#ax.set_ylim(0,500)




#print(len(df))



#print(extinctionTimes)

data08=extinctionTimes

#######################################################################################################################################
#######################################################################################################################################
varFactor = 0.9

g_x_odd = -0.05 #Initial treatment (eg ALK TKI)
g_y_odd = 0.005

b_x_odd = (g_x_odd+varFactor)/2
b_y_odd = (g_y_odd+varFactor)/2

d_x_odd = (-g_x_odd+varFactor)/2
d_y_odd = (-g_y_odd+varFactor)/2

g_x_even = 0.005 
g_y_even = -0.05 

b_x_even = (g_x_even+varFactor)/2
b_y_even = (g_y_even+varFactor)/2

           
d_x_even = (-g_x_even+varFactor)/2
d_y_even = (-g_y_even+varFactor)/2

trans1 = 0.0001 
trans2 = trans1

relativePopInit = [9000,1000];

#Check for negative rates:
for parameter in [b_x_odd,b_y_odd,d_x_odd,d_y_odd,b_y_even,b_x_even,d_x_even,d_y_even]:
    if parameter < 0:
        raise NegativeRate

ApB=relativePopInit[1]/relativePopInit[0]
        
T_initial = np.log(((g_y_odd-g_x_odd)*(g_x_even-g_x_odd)+trans1*(g_y_odd+g_x_even-g_x_odd-g_y_even))/((trans1+ApB*(trans1+g_y_odd-g_x_odd))*(g_y_odd-g_y_even)))/(trans1+g_y_odd-g_x_odd)


k=(((g_y_odd-g_x_odd)*(g_x_even-g_x_odd)+trans1*(g_y_odd+g_x_even-g_x_odd-g_y_even))*(g_y_odd-g_y_even))/(((g_y_odd-g_y_even)*(g_x_even-g_y_even)+trans1*(g_y_odd+g_x_even-g_x_odd-g_y_even))*(g_x_even-g_x_odd))

T_treat_1 = 10.0
T_treat_2 = T_treat_1/k

T_treat = T_treat_1+T_treat_2

#T_initial -> T_treat_1 -> T_treat_2 -> T_treat_1

print(T_initial,T_treat_1,T_treat_2)



    





relativePop = relativePopInit


          
#-----------------Stoch----------------------------------------

def update(state):
    
    rt,re = random(), random()
    
    timeAfter = (state[2]-T_initial)/T_treat-int((state[2]-T_initial)/T_treat)
    
    
    if timeAfter>(T_treat_1/T_treat) or (state[2]<T_initial):
        a =(b_x_odd+d_x_odd+trans1)*state[0] + (b_y_odd+d_y_odd)*state[1]  
        
        if b_x_odd*state[0]/a>re:
            state[0] += 1
        elif (d_x_odd*state[0]+b_x_odd*state[0])/a>re:
            state[0] -= 1
        elif (b_y_odd*state[1]+d_x_odd*state[0]+b_x_odd*state[0])/a>re:
            state[1] += 1
        elif (d_y_odd*state[1]+b_y_odd*state[1]+d_x_odd*state[0]+b_x_odd*state[0])/a>re:
            state[1] -= 1
        else:
            state[0] -= 1
            state[1] += 1
    else:
        a =(b_x_even+d_x_even)*state[0] + (b_y_even+d_y_even+trans2)*state[1]
        
        if b_x_even*state[0]/a>re:
            state[0] += 1
        elif (d_x_even*state[0]+b_x_even*state[0])/a>re:
            state[0] -= 1
        elif (b_y_even*state[1]+d_x_even*state[0]+b_x_even*state[0])/a>re:
            state[1] += 1
        elif (d_y_even*state[1]+b_y_even*state[1]+d_x_even*state[0]+b_x_even*state[0])/a>re:
            state[1] -= 1
        else:
            state[0] += 1
            state[1] -= 1
             
    state[2] = state[2] - (1/a)*log(rt)
        
    
    return state



replicate = 0

extinctionTimes=[]

data=[]



while (replicate < replicates):   
    
    extinct=False;
    
    state = [relativePopInit[0], relativePopInit[1], 0.0] #initial conditions for state vector(0=sensitive, 1=resistant, 2=time)    
    
    life_history = []
    life_history.append(state)
    
    while (extinct==False) & (state[2]<cutoff):
    
        #print "Time: ", state[2]," Sensitive: ", state[0]," Resistant: ", state[1]
            
        new_state = copy(state)

        state = update(new_state)    
        
        life_history.append(state)
        
        if ((state[0]==0) & (state[1]==0)):
            extinct=True
            extinctionTimes.append(state[2])
            
    
    life_history.append([0,0,cutoff])
    
    data.append(life_history)
    
    replicate += 1
    
#-----------------ODE----------------------------------------


data0=[];
data1=[];

timeStep=0.01;

for k in range(0, int(cutoff/timeStep)):
    
    if (float(k)*timeStep).is_integer():
        data0.append(relativePop[0]);
        data1.append(relativePop[1]);
        
    r01=trans1;
    r10=trans2;    
    
    
    timeAfter = (k*timeStep-T_initial)/T_treat-int((k*timeStep-T_initial)/T_treat)
    
    if timeAfter>(T_treat_1/T_treat) or (k*timeStep<T_initial):
        g0=g_x_odd;
        g1=g_y_odd;
        r01=trans1;
        r10=0.0;
    else:
        g0=g_x_even;
        g1=g_y_even;
        r01=0.0;
        r10=trans2;
    
    
        
    relativePop[0]=relativePop[0]+(relativePop[0]*g0+relativePop[1]*r10-relativePop[0]*r01)*timeStep;
    relativePop[1]=relativePop[1]+(relativePop[1]*g1+relativePop[0]*r01-relativePop[1]*r10)*timeStep;

totalODE=[]

for i,na in enumerate(data0):
    totalODE.append(data0[i]+data1[i])    
    

                   
#------------------Plots--------------------------

preDf=[]


#Bin the data into replicates:

for replicateNumber,replicateData in enumerate(data):
    
    binData=[[[] for _ in range(cutoff+1)],[[] for _ in range(cutoff+1)]]
    
    for i,repData in enumerate(replicateData):
        
        binData[0][int(replicateData[i][2])].append(int(replicateData[i][0]))
        binData[1][int(replicateData[i][2])].append(int(replicateData[i][1]))
        
    for i,binList in enumerate(binData[0]):
        
        tempData=[0,0,0]        
        
        if len(binData[0][i])==0:
            tempData[0]=binData[0][i-1][len(binData[0][i-1])-1]
            binData[0][i].append(binData[0][i-1][len(binData[0][i-1])-1])
        else:
            tempData[0]=np.mean(binData[0][i])
        if len(binData[1][i])==0:    
            tempData[1]=binData[1][i-1][len(binData[1][i-1])-1]
            binData[1][i].append(binData[1][i-1][len(binData[0][i-1])-1])
        else:
            tempData[1]=np.mean(binData[1][i])
        
        tempData[2]=sum(tempData)
        
        preDf.append([float(i),'$A_R$',str(replicateNumber),tempData[1]]) #change order here to change order of cells on tsplot
        preDf.append([float(i),'$B_R$',str(replicateNumber),tempData[0]])
        
        preDf.append([float(i),'Total',str(replicateNumber),tempData[2]])

odepreDf=[]
odepreDf2=[]

for k in range(cutoff):
    odepreDf.append([float(k),'ODE',replicates,data0[k]])
    
for k in range(cutoff):
    odepreDf2.append([float(k),'ODE2',replicates,data1[k]]) 
    
for k in range(cutoff):
    odepreDf2.append([float(k),'ODE3',replicates+1,totalODE[k]]) 

            
df=pandas.DataFrame(preDf,columns=['Time (days)', 'Cell Type','subject', 'Cell Number'])

odeDf=pandas.DataFrame(odepreDf,columns=['Time (days)', 'Cell Type','subject', 'Cell Number'])
odeDf2=pandas.DataFrame(odepreDf2,columns=['Time (days)', 'Cell Type','subject', 'Cell Number'])

#tsplot:

#sns.set_style("whitegrid", {'axes.grid' : False})

sns.set(font_scale=1.5)

ax = sns.tsplot(time="Time (days)", value="Cell Number",unit="subject", condition="Cell Type",data=df,err_style="unit_traces")

ax = sns.tsplot(time="Time (days)", value="Cell Number",unit="subject", condition="Cell Type",data=odeDf,color='k',linestyle='--',legend=True)
ax = sns.tsplot(time="Time (days)", value="Cell Number",unit="subject", condition="Cell Type",data=odeDf2,color='k',linestyle='--',legend=False)

#ax.set(yscale="log")

ax.grid(color='k',linestyle=':')

ax.set_xlim(0,200)
#ax.set_ylim(0,500)




#print(len(df))



#print(extinctionTimes)

data09=extinctionTimes

#######################################################################################################################################
#######################################################################################################################################
varFactor = 1.0

g_x_odd = -0.05 #Initial treatment (eg ALK TKI)
g_y_odd = 0.005

b_x_odd = (g_x_odd+varFactor)/2
b_y_odd = (g_y_odd+varFactor)/2

d_x_odd = (-g_x_odd+varFactor)/2
d_y_odd = (-g_y_odd+varFactor)/2

g_x_even = 0.005 
g_y_even = -0.05 

b_x_even = (g_x_even+varFactor)/2
b_y_even = (g_y_even+varFactor)/2

           
d_x_even = (-g_x_even+varFactor)/2
d_y_even = (-g_y_even+varFactor)/2

trans1 = 0.0001 
trans2 = trans1

relativePopInit = [9000,1000];

#Check for negative rates:
for parameter in [b_x_odd,b_y_odd,d_x_odd,d_y_odd,b_y_even,b_x_even,d_x_even,d_y_even]:
    if parameter < 0:
        raise NegativeRate

ApB=relativePopInit[1]/relativePopInit[0]
        
T_initial = np.log(((g_y_odd-g_x_odd)*(g_x_even-g_x_odd)+trans1*(g_y_odd+g_x_even-g_x_odd-g_y_even))/((trans1+ApB*(trans1+g_y_odd-g_x_odd))*(g_y_odd-g_y_even)))/(trans1+g_y_odd-g_x_odd)


k=(((g_y_odd-g_x_odd)*(g_x_even-g_x_odd)+trans1*(g_y_odd+g_x_even-g_x_odd-g_y_even))*(g_y_odd-g_y_even))/(((g_y_odd-g_y_even)*(g_x_even-g_y_even)+trans1*(g_y_odd+g_x_even-g_x_odd-g_y_even))*(g_x_even-g_x_odd))

T_treat_1 = 10.0
T_treat_2 = T_treat_1/k

T_treat = T_treat_1+T_treat_2

#T_initial -> T_treat_1 -> T_treat_2 -> T_treat_1

print(T_initial,T_treat_1,T_treat_2)



    





relativePop = relativePopInit


          
#-----------------Stoch----------------------------------------

def update(state):
    
    rt,re = random(), random()
    
    timeAfter = (state[2]-T_initial)/T_treat-int((state[2]-T_initial)/T_treat)
    
    
    if timeAfter>(T_treat_1/T_treat) or (state[2]<T_initial):
        a =(b_x_odd+d_x_odd+trans1)*state[0] + (b_y_odd+d_y_odd)*state[1]  
        
        if b_x_odd*state[0]/a>re:
            state[0] += 1
        elif (d_x_odd*state[0]+b_x_odd*state[0])/a>re:
            state[0] -= 1
        elif (b_y_odd*state[1]+d_x_odd*state[0]+b_x_odd*state[0])/a>re:
            state[1] += 1
        elif (d_y_odd*state[1]+b_y_odd*state[1]+d_x_odd*state[0]+b_x_odd*state[0])/a>re:
            state[1] -= 1
        else:
            state[0] -= 1
            state[1] += 1
    else:
        a =(b_x_even+d_x_even)*state[0] + (b_y_even+d_y_even+trans2)*state[1]
        
        if b_x_even*state[0]/a>re:
            state[0] += 1
        elif (d_x_even*state[0]+b_x_even*state[0])/a>re:
            state[0] -= 1
        elif (b_y_even*state[1]+d_x_even*state[0]+b_x_even*state[0])/a>re:
            state[1] += 1
        elif (d_y_even*state[1]+b_y_even*state[1]+d_x_even*state[0]+b_x_even*state[0])/a>re:
            state[1] -= 1
        else:
            state[0] += 1
            state[1] -= 1
             
    state[2] = state[2] - (1/a)*log(rt)
        
    
    return state



replicate = 0

extinctionTimes=[]

data=[]



while (replicate < replicates):   
    
    extinct=False;
    
    state = [relativePopInit[0], relativePopInit[1], 0.0] #initial conditions for state vector(0=sensitive, 1=resistant, 2=time)    
    
    life_history = []
    life_history.append(state)
    
    while (extinct==False) & (state[2]<cutoff):
    
        #print "Time: ", state[2]," Sensitive: ", state[0]," Resistant: ", state[1]
            
        new_state = copy(state)

        state = update(new_state)    
        
        life_history.append(state)
        
        if ((state[0]==0) & (state[1]==0)):
            extinct=True
            extinctionTimes.append(state[2])
            
    
    life_history.append([0,0,cutoff])
    
    data.append(life_history)
    
    replicate += 1
    
#-----------------ODE----------------------------------------


data0=[];
data1=[];

timeStep=0.01;

for k in range(0, int(cutoff/timeStep)):
    
    if (float(k)*timeStep).is_integer():
        data0.append(relativePop[0]);
        data1.append(relativePop[1]);
        
    r01=trans1;
    r10=trans2;    
    
    
    timeAfter = (k*timeStep-T_initial)/T_treat-int((k*timeStep-T_initial)/T_treat)
    
    if timeAfter>(T_treat_1/T_treat) or (k*timeStep<T_initial):
        g0=g_x_odd;
        g1=g_y_odd;
        r01=trans1;
        r10=0.0;
    else:
        g0=g_x_even;
        g1=g_y_even;
        r01=0.0;
        r10=trans2;
    
    
        
    relativePop[0]=relativePop[0]+(relativePop[0]*g0+relativePop[1]*r10-relativePop[0]*r01)*timeStep;
    relativePop[1]=relativePop[1]+(relativePop[1]*g1+relativePop[0]*r01-relativePop[1]*r10)*timeStep;

totalODE=[]

for i,na in enumerate(data0):
    totalODE.append(data0[i]+data1[i])    
    

                   
#------------------Plots--------------------------

preDf=[]


#Bin the data into replicates:

for replicateNumber,replicateData in enumerate(data):
    
    binData=[[[] for _ in range(cutoff+1)],[[] for _ in range(cutoff+1)]]
    
    for i,repData in enumerate(replicateData):
        
        binData[0][int(replicateData[i][2])].append(int(replicateData[i][0]))
        binData[1][int(replicateData[i][2])].append(int(replicateData[i][1]))
        
    for i,binList in enumerate(binData[0]):
        
        tempData=[0,0,0]        
        
        if len(binData[0][i])==0:
            tempData[0]=binData[0][i-1][len(binData[0][i-1])-1]
            binData[0][i].append(binData[0][i-1][len(binData[0][i-1])-1])
        else:
            tempData[0]=np.mean(binData[0][i])
        if len(binData[1][i])==0:    
            tempData[1]=binData[1][i-1][len(binData[1][i-1])-1]
            binData[1][i].append(binData[1][i-1][len(binData[0][i-1])-1])
        else:
            tempData[1]=np.mean(binData[1][i])
        
        tempData[2]=sum(tempData)
        
        preDf.append([float(i),'$A_R$',str(replicateNumber),tempData[1]]) #change order here to change order of cells on tsplot
        preDf.append([float(i),'$B_R$',str(replicateNumber),tempData[0]])
        
        preDf.append([float(i),'Total',str(replicateNumber),tempData[2]])

odepreDf=[]
odepreDf2=[]

for k in range(cutoff):
    odepreDf.append([float(k),'ODE',replicates,data0[k]])
    
for k in range(cutoff):
    odepreDf2.append([float(k),'ODE2',replicates,data1[k]]) 
    
for k in range(cutoff):
    odepreDf2.append([float(k),'ODE3',replicates+1,totalODE[k]]) 

            
df=pandas.DataFrame(preDf,columns=['Time (days)', 'Cell Type','subject', 'Cell Number'])

odeDf=pandas.DataFrame(odepreDf,columns=['Time (days)', 'Cell Type','subject', 'Cell Number'])
odeDf2=pandas.DataFrame(odepreDf2,columns=['Time (days)', 'Cell Type','subject', 'Cell Number'])

#tsplot:

#sns.set_style("whitegrid", {'axes.grid' : False})

sns.set(font_scale=1.5)

ax = sns.tsplot(time="Time (days)", value="Cell Number",unit="subject", condition="Cell Type",data=df,err_style="unit_traces")

ax = sns.tsplot(time="Time (days)", value="Cell Number",unit="subject", condition="Cell Type",data=odeDf,color='k',linestyle='--',legend=True)
ax = sns.tsplot(time="Time (days)", value="Cell Number",unit="subject", condition="Cell Type",data=odeDf2,color='k',linestyle='--',legend=False)

#ax.set(yscale="log")

ax.grid(color='k',linestyle=':')

ax.set_xlim(0,200)
#ax.set_ylim(0,500)




#print(len(df))



#print(extinctionTimes)

data10=extinctionTimes

#######################################################################################################################################
#######################################################################################################################################