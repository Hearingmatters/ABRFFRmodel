close all
clear all

Fs=100000; % if you want to change this value for the model, then change this hardcoded in the cochlear_model.py

f0 = 1000;
ramp = 1;
ton=10e-3;
tsig=250e-3;
tsin=0:1/Fs:tsig-(1/Fs);
tend=250e-3;

Wind=MakeWindow(numel(tsin),ton,Fs,1);
sinst=sin(2*pi*f0*tsin).*Wind;

if(ramp)
    tone =[0.999*sinst zeros(1,round(tend*Fs)-numel(sinst))]; 
    %rampsoundhanning(sin(2*pi*f0*t), Fs, rise);
else
    tone = sin(2*pi*f0*t);
end

spl = 60;

channels = numel(spl);
%for k = 1:channels
%    stim(k,:) = tone;
%end

for k = 1:channels
    stim(k,:) = sqrt(2)*tone./(sqrt(mean(tone.^2)));
    stim=[stim  zeros(1,round(200e-3*Fs))  stim zeros(1,round(100e-3*Fs))];
end

normalizeRMS=zeros(1,channels);
%the model is now hardcoded to multiply the input by 20e-6*10(spl/20)
%if you want your spl to correspond to rms, you need to do the above
%compensation..
%similarly when running a condensation click 0-1, and you need peSPL, you
%need to multipy the 1 of your click by sqrt(2)*2 to have a ptp value
%corresponding to a pure-tone of the same rms

subject=1;
irregularities=ones(1,channels);
%to generate reflection-source emissions (0 or 1)

sheraPo=0.0610; 
%This value is not used here, and overwritten by the values of the 'StartingPoles.dat' files in the sysfiles folder
%If you want to have a sloping hearing loss rename "StartingPoles.dat" to
%"StartingPolesNH" and "rename StartingPolesHI.dat" to "StartingPoles".

probes = 'all';
clear tone t rise dur f0
save input;


