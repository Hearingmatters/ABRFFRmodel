function CosWindow=MakeWindow(sigdur,onset,fs,both)

%onset=10e-3;
%fs=48800;
%sigdur=6000; %in samples

Wdur=2*round(onset*fs);
W=hanning(Wdur);

if both==0
    CosWindow=[W(1:Wdur/2)' ones(1,sigdur-Wdur/2,1)];
else
    CosWindow=[W(1:Wdur/2)' ones(1,sigdur-Wdur,1) W(Wdur/2+1:end)'];
end