import numpy as np
import scipy as sp
from scipy import signal
import scipy.io as sio
import cochlear_model
import os
#def erb2hz(x):
#    return 229*(10.**(x/21.3)-1)
#def next_pow2(x):
#    return 2**(np.ceil(np.log2(x)))
import multiprocessing as mp
import ctypes as c
import time

#import sys
#print (sys.version)
Oversampling=1
sectionsNo=1000
p0=float(2e-5)
par=sio.loadmat('/home/gmehraei/ABB_model/StimInput/inputclick.mat')

#probes=np.array(par['probes']) 

probes = np.array([0,9])  # what conditions to run
probe_points = probes
Fs=par['Fs']
Fs=Fs[0][0]
stim=par['stim']
stim = np.float64(stim)  # the input is somtimes read uint8 so convert
stim=stim[probes]
# print stim
normalizeRMS=par['normalizeRMS']
#normalizeRMS=np.array(normalizeRMS)
#normalizeRMS=normalizeRMS[probes]
spl=par['spl'] 
spl=np.array(spl[0])
spl=spl[probes]
print spl

channels=np.size(spl)
#channels=par['channels']
#channels=channels[1][0]
subjectNo=int(par['subject'])
lgt=len(stim[0])
norm_factor=p0*10.**(spl/20.)

#sheraPo=0.06
sheraPo = np.loadtxt('/home/gmehraei/ABRFFRmodel/sysfiles/StartingPoles.dat', delimiter=',')
sheraPo = np.array(sheraPo)
irr_on=np.array(par['irregularities'])
d=len(stim[0].transpose())
print("running cochlear simulation")
print spl

def solve_one_cochlea(model): #definition here, to have all the parameter implicit
    #i=model[3]
    coch=model[0]
    coch.init_model(model[1],Oversampling*Fs,sectionsNo,probe_points,Zweig_irregularities=model[2],sheraPo=sheraPo,subject=subjectNo) #model needs to be init here because if not pool.map crash
    coch.solve()
#    f=open("out/v"+str(i+1)+".np",'wb')
#    np.array(coch.Vsolution,dtype='=d').tofile(f)
#    f.close()
#    f=open("out/y"+str(i+1)+".np",'wb')
#    np.array(coch.Ysolution,dtype='=d').tofile(f)
#    f.close()
#    f=open("out/E"+str(i+1)+".np",'wb')
#    np.array(coch.oto_emission,dtype='=d').tofile(f)
#    f.close()
#    f=open("out/F"+str(i+1)+".np",'wb')
#    np.array(coch.cf,dtype='=d').tofile(f)
#    f=open("out/Fc"+str(i+1)+".np",'wb')
#    np.array(coch.f_resonance,dtype='=d').tofile(f)
#    f.close()
    return [coch.Vsolution, coch.Ysolution, coch.oto_emission,
            coch.stim[0:len(coch.oto_emission)], coch.f_resonance]

for i in range(channels):
    #for clicks
    stimRms=1/(np.sqrt(2));
    #stimRms=1
    stim[i]=stim[i]/stimRms*norm_factor[i]
    
sig=stim

cochlear_list=[ [cochlear_model.cochlea_model(),sig[i],irr_on[0][i],i] for i in range(channels)]


if __name__ == "__main__":
    p=mp.Pool(int(mp.cpu_count()),maxtasksperchild=1)
    result=p.map(solve_one_cochlea,cochlear_list)
   
    Vresult = np.ndarray([len(result[0][0].transpose()), sectionsNo + 1, channels])
    Yresult = np.ndarray([len(result[0][0].transpose()), sectionsNo + 1, channels])
    Emission = np.zeros([len(result[0][0].transpose()), channels])
    Resampled_stimulus = np.zeros([len(result[0][3].transpose()), channels])
    Fc = result[0][4]
    for i in range(channels):
        Vresult[:, :, i] = result[i][0].transpose()
        Yresult[:, :, i] = result[i][1].transpose()
        Emission[:, i] = result[i][2]
        Resampled_stimulus[:, i] = result[i][3]


    sio.savemat("/home/gmehraei/ABB_model/CochleaOutput/outputclick_test.mat",
             mdict={"Velocity": Vresult, "Displacement": Yresult,
                   "OtoAcousticEmission": Emission,
                   "OutStimulus": Resampled_stimulus,
                   "model_sample_rate": float(Fs * Oversampling), "Fc": Fc})

    
    p.close()
    p.join()

    
