close all;
clear all;
%*************************************************************************
Fs=100000;
type='NH';

name=[type,'ClicksMEL'];
%move the correct poles to the execution folder
if type=='HI'
    copyfile('../sysfiles/HIPoles/StartingPoles.dat','../sysfiles')
elseif type=='NH'
    copyfile('../sysfiles/NHPoles/StartingPoles.dat','../sysfiles')
end

% Model params
channels = 11;    %stimulus blocks
normalizeRMS=zeros(1,channels);
subject=1;
spl =[0:10:100];
irregularities=ones(1,channels);
sheraPo=0.0610; 
probes = 'all';

%% click generation
for k=1:channels
    Cdur=round(80e-6*Fs);
    SC=[zeros(1,round(Fs*1e-3))  ones(1,Cdur) zeros(1,round(Fs*50e-3))];
    stim(k,:)=2*sqrt(2)*SC; 
end

save('input.mat','stim','Fs','spl','channels','normalizeRMS','subject','irregularities','sheraPo','probes')
system('python run_cochlear_modelold.py')   
%movefile('output.mat',['./out/Clicks/',name,'.mat'])

break
