close all
clear all

%Loads in the BM velocity simulations from the model
%path='/media/tic/ModellingTLM/out/ClickNH_rms/'; %should end with a /
%path='/media/tic/ModellingTLM/out/ClickLCHI_oldVth/'; %should end with a /
%file='probingClickSCL'; %should be the file name without the variable in it
%L=0:10:100;
L=1;
%load Fcnew.mat %the probing frequencies
load outputtest_new.mat
FS=100000;
implnt=0;

%% parameters

% for the auditory nerve part
nrep=1; %number of stimulus repetitions

%for Vbm to Ycilia stage: with vbm this is a LP filter, for ybm
%this would be a HP filter
Mvel=41e-6; %from 1 kHz max BM vel from Heinz Stimuli 100dB.. %0.09645;, where the fuck does this high value come from.. %Max BM velocity in the model
Mu=200e-9;  %Max cilia displacement in the model
Fcut=513;   %this value should come somewhere from Sellick et al.
tau_c=1/(2*pi*Fcut);
Fgain=Mu/Mvel;
%FA=Fgain/(1+tau_c*2*FS);
%FB=(tau_c*2*FS-1)/(tau_c*2*FS+1);

%for the S shape IHC NL model that fits to the data of Fig12 in LopezPoveda
%2006 using the Russel and Sellick 1983 300Hz IHC receptor potential I/O.
%Note that axis has changed from stim pressure to be -200 to +200nm (i.e., hair bundle displacement)
Amp=1.28e-3;
Off=-0.26e-3;
beta=0.075/1e-9; %slope of the function
alpha=(1/beta)*log((Amp+Off)/-Off);


%for the IHC low-pass filter a la Zilany et al., with different
%cut-off here 2.5 iso 3800/3000kHz in Zilany!
F_LPC=1000; % 2500; %4500;     %4500 in Heinz model! %3800/3000 in zilani in original  %IHC lowpass cutoff frequency, Heinz uses 4500Hz 
LPk=2 ; %7;         %order of the lowpass filter
C=2*FS;
C1LP=(C-(2*pi*F_LPC))/(C+(2*pi*F_LPC));
C2LP=(2*pi*F_LPC)/((2*pi*F_LPC)+C);
gain=1;
sec=1000;

%% Load in the stimulus here:
for rL=2
    %for NS=1:26 %for all cochlear sections (26 files to simulate the whole BM)
       Vbm_tmp=squeeze(Velocity(:,2:end,rL));
        %% load in Vbm and resample
        %display([num2str(NS),'/26'])
        %eval(['load ',path,file,num2str(L(rL)),'_',num2str(NS),'.dat'])
    
        %determine the storing pointer, because there were some sections
        %missing in the simulations (don't worry about this, it is fine)
        %if NS<14; %19 probes
         %   sec=19;
         %   ns=(NS-1)*sec; %pointer-1 of where to start writing next
        %elseif NS>=14 && NS<26 %20 probes
        %    sec=20; %number of sections simulated in one probing simulation
        %    ns=(13*19)+(NS-14)*sec;
        %elseif NS==26
         %   sec=11; %number of sections simulated in one probing simulation
          %  ns=13*19+12*20; 
        %end

        %clear s
        %for n=1:sec
        %    eval(['s(:,',num2str(n),')=',file,num2str(L(rL)),'_',num2str(NS),'(:,',num2str(sec+n),');']) %use BMvel!       
        %end
        ns=0;
        Vbm=resample(Vbm_tmp,FS,model_sample_rate); %resampling from old to new FS
        %eval(['clear ',file,num2str(L(rL)),'_',num2str(NS)]) %get rid of the originally loaded vector   
        %Vbm=Vbm.*1000; %temp. 
        %% do calculations for each simulated section
        for n=1:sec                
            for k=1:size(Vbm,1);
%                 %% Vbm to YCilia filter a la Shamma 1986
%                 if k==1;
%                     yc(k)=0;
%                 else
%                     yc(k)=FA*Vbm(k,n)+FA*Vbm(k-1,n)+FB*yc(k-1);
%                 end
                if Fc(n)<Fcut;
                    yc(k,n)=Fgain*Vbm(k,n); %low freqs proportional to Vbm
                else
                    yc(k,n)=Fgain*Vbm(k,n); %high freqs proporional to Ybm 
                    %FA*Vbm(k,n)+FA*Vbm(k-1,n)+FB*yc(k-1,n);
                end
            
                %% Ycilia to receptor potential a la Mountain & Hubbard 1996 book to match Russel et al. 1986 data
                %Popen(k)=1./((1+exp(-(yc(k)-x0)/S0))*(1+exp(-(yc(k)-x1)/S1)));
                %G0=0; %just a test here..with 0
                %VihcNF(k)=Gmax*Popen(k)+G0; %non-filtered IHC receptor potential
                VihcNF(k)=Off+Amp*(1./(1+exp(beta*(alpha-yc(k)))));
            end
            
            IHC1=0*ones(LPk+1,1); 
            IHC2=0*ones(LPk+1,1); 
            for k=1:size(Vbm,1);
                %% IHC Low-pass filter
                IHC1(1)=gain*VihcNF(k);
                    for r=1:LPk
                        IHC1(r+1)=C1LP*IHC2(r+1)+C2LP*(IHC1(r)+IHC2(r));
                    end

                    for r=1:(LPk+1)
                        IHC2(r)=IHC1(r);
                    end
                Vihc(k)=IHC1(LPk+1); 
                
            end %1 section over time
                
            %innerHC=repmat(Vihc,1,nrep); %do the repetitions here as an input to the AN model
            %% call the auditory nerve model
             CF=Fc(n); %Char Freq of the synapse
             if CF>80; %the AN model only works for freq higher than 80 Hz
                fiberType = 1; %Low spont
                [synoutLS,psthLS] = Verhulst2014_NOFD(Vihc,CF,nrep,1/FS,fiberType,implnt);
                fiberType = 2; %Med spont
                [synoutMS,psthMS] = Verhulst2014_NOFD(Vihc,CF,nrep,1/FS,fiberType,implnt);
                fiberType = 3; %High spont 
                [synoutHS,psthHS] = Verhulst2014_NOFD(Vihc,CF,nrep,1/FS,fiberType,implnt);
             else 
                synoutLS=NaN(size(Vihc,1),1);
                synoutMS=NaN(size(Vihc,1),1);
                synoutHS=NaN(size(Vihc,1),1);
                psthLS=NaN(size(Vihc,1),1);
                psthMS=NaN(size(Vihc,1),1);
                psthHS=NaN(size(Vihc,1),1);
             end            
            
            % save the responses
            BM(:,ns+n)=Vbm(:,n);
            YC(:,ns+n)=yc(:,n);
            VIHCNF(:,ns+n)=VihcNF;
            VIHC(:,ns+n)=Vihc;
            ANLS(:,ns+n)=synoutLS;
            ANMS(:,ns+n)=synoutMS;
            ANHS(:,ns+n)=synoutHS;      
               
        end %end of all sections
    end %end of all NS
    
    %save the whole cochlear partition
    eval(['save(''./out/ClicksBBANHI/BM',num2str(L(rL)),'.mat'',''BM'')'])
    eval(['save(''./out/ClicksBBANHI/YC',num2str(L(rL)),'.mat'',''YC'')'])
    eval(['save(''./out/ClicksBBANHI/VIHCNF',num2str(L(rL)),'.mat'',''VIHCNF'')'])
    eval(['save(''./out/ClicksBBANHI/VIHC',num2str(L(rL)),'.mat'',''VIHC'')'])
    eval(['save(''./out/ClicksBBANHI/ANLS',num2str(L(rL)),'.mat'',''ANLS'')'])
    eval(['save(''./out/ClicksBBANHI/ANMS',num2str(L(rL)),'.mat'',''ANMS'')'])
    eval(['save(''./out/ClicksBBANHI/ANHS',num2str(L(rL)),'.mat'',''ANHS'')'])
%end %end of all different levels

