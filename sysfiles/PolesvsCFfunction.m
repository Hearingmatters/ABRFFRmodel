close all
clear all

Nsec=1000;
BMlength=35e-3 - 1e-3; %35 mm BM length - 1 mm helicotrema width
dx=BMlength/Nsec;
Ga=20682;
Galpha=61.765;
Gbeta=140.6;
%Greenwoodmap
x=dx:dx:BMlength;
fres=(Ga*10.^(-Galpha*x))-Gbeta; %the greenwood function

%%parameters: comes from fit to tuning data
acat=11.4596;
bcat=0.2533;

aQ=0.3570;
bQ=-0.7933;

a=0.2927;
b=-1.2337;

Qhuman=acat*(fres/1000).^bcat;
%calculation of the starting poles using the powerlaw fit of the Qerb/pole
%function of Verhulst et al. 2012
StartingPoles=aQ*Qhuman.^bQ;
%for stability and memory purposes in the model, do not go beyond the range
%of 0.3 and 0.03 for the model poles. So limit these
Ns=find(StartingPoles<0.037); %originally 0.037
StartingPoles(Ns)=0.037; %originally 0.037

%uncomment this last bit if you want to write a new startingpoles file
%fid = fopen('StartingPoles.dat','w');
%fprintf(fid,'%E\n',StartingPoles);
%fclose(fid);

figure,plot(StartingPoles)
xlabel('sectionNo'),ylabel('PoleLocation')

figure,plot(fres/1000,StartingPoles)
xlabel('CF'),ylabel('PoleLocation')

%model Q calculated from the StartingPoles
ModelQ=a*StartingPoles.^b;
figure,plot(fres/1000,ModelQ)
xlabel('CF'),ylabel('Qerb')

