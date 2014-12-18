%clear all;
%mex -v zilany2009_humanized1_IHC.c complex.c  
clear all;
%mex -f /mnt/localhd/home/save/Documents/MATLAB/mexopts.sh

%OL PC
%mex -f /home/sarah/Documents/MATLAB/mexopts.sh -v zilany2009_humanized1_SynapseHeinz.c complex.c 
%mex -f /home/sarah/Documents/MATLAB/mexopts.sh -v zilany2009_NOFD.c complex.c
%mex -f /home/sarah/Documents/MATLAB/mexopts.sh -v zilany2009_NOFD_PLA.c complex.c
%mex -f /home/sarah/Documents/MATLAB/mexopts.sh -v zilany2009_NOFD_noMatch.c complex.c

%mex -f /home/staralfur/Documents/MATLAB/mexopts.sh -v Verhulst2014_NOFD.c complex.c
mex -f /home/staralfur/Documents/MATLAB/mexopts.sh -v Verhulst2014_NOFD_TH.c complex.c
mex -f /home/staralfur/Documents/MATLAB/mexopts.sh -v model_Synapse.c complex.c
%mex -f /home/staralfur/Documents/MATLAB/mexopts.sh -v model_Synapse_CI.c complex.c
%mex -f /home/sarah/Documents/MATLAB/mexopts.sh -v Verhulst2014_NOFD_PLA.c complex.c
%mex -f /home/staralfur/Documents/MATLAB/mexopts.sh -v Verhulst2014_NOFD_PLA_ffGN.c complex.c
%mex -f /home/sarah/Documents/MATLAB/mexopts.sh -v Verhulst2014_NOFD_Goldie.c complex.c

%mex -f /home/staralfur/Documents/MATLAB/mexopts.sh -v Verhulst2014.c complex.c
%mex -f /home/staralfur/Documents/MATLAB/mexopts.sh -v Verhulst2014_PLA.c complex.c
%mex -f /home/staralfur/Documents/MATLAB/mexopts.sh -v Verhulst2014_PLA_ffGN.c complex.c

%mex -f /home/staralfur/Documents/MATLAB/mexopts.sh -v AltoeAN.c complex.c
