%function ANF(name)
close all
clear all

L=0:10:100;
%name='NHClicksME';
name='output';

FS=100000;
implnt=0;
n=1;

tic
load(['../out/Clicks/',name,'.mat'],'Velocity','Fc');
CF=Fc(2:2:numel(Fc));
%% parameters
nrep=1; %number of stimulus repetitions
Mvel=41e-6; %from 1 kHz max BM vel from Heinz Stimuli 100dB.. %0.09645;, where the fuck does this high value come from.. %Max BM velocity in the model
Mu=200e-9;  %Max cilia displacement in the model
Fcut=513;   %this value should come somewhere from Sellick et al.
tau_c=1/(2*pi*Fcut);
Fgain=Mu/Mvel;

%for the IHC low-pass filter a la Zilany et al., with different
%cut-off here 2.5 iso 3800/3000kHz in Zilany!
F_LPC=1000; % 2500; %4500;     %4500 in Heinz model! %3800/3000 in zilani in original  %IHC lowpass cutoff frequency, Heinz uses 4500Hz
LPk=2 ; %7;         %order of the lowpass filter
C=2*FS;
C1LP=(C-(2*pi*F_LPC))/(C+(2*pi*F_LPC));
C2LP=(2*pi*F_LPC)/((2*pi*F_LPC)+C);
gain=1;

%% do for each stimulus level
for m=9
    display(num2str(m))
    %% do calculations for each simulated section
    for n=2:2:numel(Fc) %do for every other section
        display(num2str(n/2))
        %% IHC deflection and nonlinearity
        for k=1:size(Velocity,1);
            yc(k)=Fgain*Velocity(k,n,m);
            %VihcNF(k)=Off+Amp*(1./(1+exp(beta*(alpha-yc(k)))));
            %try the old nonlinearity
            A0=0.0008;       %0.1 scalar in IHC nonlinear function
            B=2000*6000;  %2000 par in IHC nonlinear function
            C=0.33;             %1.74 par in IHC nonlinear function
            D=200e-9;         %6.87e-9; %par in IHC nonlinear function
            if yc(k)>=0
                Apos=A0;
                VihcNF(k)=Apos.*log(1+B*abs(yc(k)));
            else
                Aneg=-A0*(((abs(yc(k)).^C)+D)./((3*abs(yc(k)).^C)+D));
                VihcNF(k)=Aneg.*log(1+B*abs(yc(k)));
            end
        end
        
        %% IHC Low-pass filter
        IHC1=0*ones(LPk+1,1);
        IHC2=0*ones(LPk+1,1);
        for k=1:size(Velocity,1);
            IHC1(1)=gain*VihcNF(k);
            for r=1:LPk
                IHC1(r+1)=C1LP*IHC2(r+1)+C2LP*(IHC1(r)+IHC2(r));
            end
            
            for r=1:(LPk+1)
                IHC2(r)=IHC1(r);
            end
            Vihc(k)=IHC1(LPk+1);
        end %1 section over time
        
        %% call the auditory nerve model
        if Fc(n)>80; %the AN model only works for freq higher than 80 Hz
            fiberType = 1; %Low spont
            [ANLS,psthLS] = Verhulst2014_NOFD_TH(Vihc,Fc(n),nrep,1/FS,fiberType,implnt);
            fiberType = 2; %Med spont
            [ANMS,psthMS] = Verhulst2014_NOFD_TH(Vihc,Fc(n),nrep,1/FS,fiberType,implnt);
            fiberType = 3; %High spont
            [ANHS,psthHS] = Verhulst2014_NOFD_TH(Vihc,Fc(n),nrep,1/FS,fiberType,implnt);
        else
            ANLS=NaN(size(Vihc,1),1);
            ANMS=NaN(size(Vihc,1),1);
            ANHS=NaN(size(Vihc,1),1);
            psthLS=NaN(size(Vihc,1),1);
            psthMS=NaN(size(Vihc,1),1);
            psthHS=NaN(size(Vihc,1),1);
        end
        IHC(:,n/2)=Vihc;
        LS(:,n/2)=ANLS;
        MS(:,n/2)=ANMS;
        HS(:,n/2)=ANHS;
    end %end of all CFs
    cd '/home/gmehraei/ABRFFRmodel/out/Clicks/'
    
    save(['TH_IHC_',num2str(L(m)),'.mat'],'IHC')
    save(['TH_ANLS_',num2str(L(m)),'.mat'],'LS')
    save(['TH_AMLS_',num2str(L(m)),'.mat'],'MS')
    save(['TH_ANHS_',num2str(L(m)),'.mat'],'HS')
    
    
end %end for all levels
%matlabpool(close)

