import numpy as np
from numpy import fft
import scipy as sp
from scipy import signal
import scipy.io as sio
from scipy.io import wavfile
import matplotlib.pyplot as plt
import main_loop_update
from os import sys
def erb2hz(x):
    return 229*(10.**(x/21.3)-1)
def next_pow2(x):
    return 2**(np.ceil(np.log2(x)))
import multiprocessing
import time

Oversampling=1
sectionsNo=1000
p0=float(2e-5)
par=sio.loadmat('input1kv2.mat')


probes=np.array([1,100]) #Hack, fix me!
probe_points=probes
Fs=par['Fs']
Fs=Fs[0][0]
stim=par['stim']
stim=np.float64(stim) #the input is somtimes read uint8 so convert 
#print stim
spl=par['spl'] #the reference is the left channel
spl=np.array(spl[0])
channels=par['channels']

channels=channels[0][0]
subjectNo=int(par['subject'])
lgt=len(stim[0])
norm_factor=p0*10.**(spl/20.)

print norm_factor

sheraPo=np.loadtxt('StartingPoles.dat',delimiter=',')
sheraPo=np.array(sheraPo)
#print sheraPo 
#sheraPo=np.append([0.037],sheraPo)

#sheraPo=par['sheraPo']
#sheraPo=np.array(sheraPo[0])
irr_on=np.array(par['irregularities'])
print(irr_on)
def solve_one_cochlea(model): #definition here, to have all the parameter implicit
   
    model[0].init_model(model[1],Oversampling*Fs,sectionsNo,probe_points,Zweig_irregularities=model[2],sheraPo=sheraPo,subject=subjectNo) #model needs to be init here because if not map crash
    model[0].solve()
    return [model[0].Vsolution,model[0].Ysolution,model[0].oto_emission,model[0].stim[0:len(model[0].oto_emission)],model[0].f_resonance]

for i in range(channels):
    #stimRms=np.sqrt(sum(stim[i]**2)/len(stim[i]))
    #stimRms=1/(2*np.sqrt(2)) #rms of click 
    stimRms=1
    stim[i]=stim[i]/stimRms*norm_factor[i]
    #print stim[2]
sig=stim
#print (sig[i])

cochlear_list=[ [main_loop_update.cochlea_model(),sig[i],irr_on[0][i]] for i in range(channels)]

print(sheraPo)
print(spl)
print(norm_factor)

#return

p=multiprocessing.Pool(int(channels))

result=p.map(solve_one_cochlea,cochlear_list)
Vresult=np.ndarray([len(result[0][0].transpose()),sectionsNo+1,channels])
Yresult=np.ndarray([len(result[0][0].transpose()),sectionsNo+1,channels])
Emission=np.zeros([len(result[0][0].transpose()), channels])
Resampled_stimulus=np.zeros([len(result[0][3].transpose()), channels])
Fc=result[0][4]
for i in range(channels):
    Vresult[:,:,i]=result[i][0].transpose()
    Yresult[:,:,i]=result[i][1].transpose()
    Emission[:,i]=result[i][2]
    Resampled_stimulus[:,i]=result[i][3]


    
sio.savemat("outputPT1_update.mat",mdict={"Velocity":Vresult,"Displacement":Yresult,"OtoAcousticEmission":Emission,"OutStimulus":Resampled_stimulus,"model_sample_rate":float(Fs*Oversampling),"Fc":Fc})

    
