import numpy as np
import scipy.io as sio
import cochlear_model
import multiprocessing as mp

Oversampling = 1
sectionsNo = 1000
p0 = float(2e-5)

# Input parameters are loaded from a mat file
par = sio.loadmat('input.mat')

probes = np.array(par['probes'])
probe_points = probes
Fs = par['Fs']
Fs = Fs[0][0]
stim = par['stim']
stim = np.float64(stim)  # the input is somtimes read uint8 so convert
spl = par['spl']
spl = np.array(spl[0])
channels = par['channels']
channels = channels[0][0]
subjectNo = int(par['subject'])
lgt = len(stim[0])
norm_factor = p0 * 10. ** (spl / 20.)

# sheraPo=0.06
sheraPo = np.loadtxt('../sysfiles/StartingPoles.dat', delimiter=',')
sheraPo = np.array(sheraPo)
irr_on = np.array(par['irregularities'])
d = len(stim[0].transpose())
print("running cochlear simulation")


#definition here, to have all the parameter implicit
def solve_one_cochlea(model):
    # i=model[3]
    coch = model[0]
    #model needs to be init here because if not pool.map crash
    coch.init_model(model[1], Oversampling * Fs, sectionsNo, probe_points,
                    Zweig_irregularities=model[2], sheraPo=sheraPo,
                    subject=subjectNo)
    coch.solve()
    return [coch.Vsolution, coch.Ysolution, coch.oto_emission,
            coch.stim[0:len(coch.oto_emission)], coch.f_resonance]

for i in range(channels):
    # stimRms=1/(2*np.sqrt(2));
    stimRms = 1
    stim[i] = stim[i] / stimRms * norm_factor[i]

sig = stim

cochlear_list = [[cochlear_model.cochlea_model(), sig[i], irr_on[0][i], i]
                 for i in range(channels)]


if __name__ == "__main__":
    p = mp.Pool(int(mp.cpu_count()), maxtasksperchild=1)
    result = p.map(solve_one_cochlea, cochlear_list)

    Vresult = np.ndarray(
        [len(result[0][0].transpose()), sectionsNo + 1, channels])
    Yresult = np.ndarray(
        [len(result[0][0].transpose()), sectionsNo + 1, channels])
    Emission = np.zeros([len(result[0][0].transpose()), channels])
    Resampled_stimulus = np.zeros([len(result[0][3].transpose()), channels])
    Fc = result[0][4]
    for i in range(channels):
        Vresult[:, :, i] = result[i][0].transpose()
        Yresult[:, :, i] = result[i][1].transpose()
        Emission[:, i] = result[i][2]
        Resampled_stimulus[:, i] = result[i][3]

    sio.savemat('output.mat',
                mdict={"Velocity": Vresult, "Displacement": Yresult,
                       "OtoAcousticEmission": Emission,
                       "OutStimulus": Resampled_stimulus,
                       "model_sample_rate": float(Fs * Oversampling),
                       "Fc": Fc})

    p.close()
    p.join()
