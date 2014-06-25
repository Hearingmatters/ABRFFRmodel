%Spike generator & Hist

close all
clear all

FS=100000;
nrep=100; %AN repetitions 120 ms pure-tone chunk
binw=0.1e-3;
WhichF=4;
Treps=1;

rng(103) %set the seed for the random number generator, so that it is always the same sequence.
jit=round(rand(1,nrep)*10e-3*FS);
%calculate the onset times of the epochs (jittering was included for AN stage)
 inner = zeros(1,200000);
 inner(1) = 1;
        onsets=inner(1,1:11000+jit(1));
        for m=2:nrep
            onsets=[onsets inner(1,1:11000+jit(m))]; %jittered stimulus in clicktrain
            %corresponds to 100kHz FS (as original ANHS)
        end        
onsettimesSamps = find(onsets)/(FS); %here we resample the onsets to match 20kHz and make in time

for L=0:10:100;
    %tic   
        eval(['load ANHS',num2str(WhichF),num2str(L),'.mat'])
        eval(['load ANLS',num2str(WhichF),num2str(L),'.mat'])

     for ch=1:5
        display(['channel no: ',num2str(ch),''])
        spiketimesHS=[];
        spiketimesLS=[];
        for n=1:Treps
            spiketimesHS =[spiketimesHS getSpikesFromRates(ANHS(:,ch),2,FS)]; %2 for HS
            spiketimesLS =[spiketimesLS getSpikesFromRates(ANLS(:,ch),1,FS)]; %1 for LS
            display(['repno: ',num2str(n),'/10'])
        end

        t=[0:1:length(ANHS)-1]/FS;

        for n=1:nrep-1
            indHS=find(spiketimesHS>=onsettimesSamps(n) & spiketimesHS<onsettimesSamps(n+1)); %here find times that fall within the nth epoch
            spiketimesHS(indHS)=spiketimesHS(indHS)-onsettimesSamps(n); %replace absolute times to times relative to onset of the epoch (so that each epoch starts at 0)       
            indLS=find(spiketimesLS>=onsettimesSamps(n) & spiketimesLS<onsettimesSamps(n+1)); %here find times that fall within the nth epoch
            spiketimesLS(indLS)=spiketimesLS(indLS)-onsettimesSamps(n); %replace absolute times to times relative to onset of the epoch (so that each epoch starts at 0)       
        end

        psthHS=histc(spiketimesHS,0:binw:500e-3)/(binw*nrep*Treps);
        psthLS=histc(spiketimesLS,0:binw:500e-3)/(binw*nrep*Treps);

        tpsth=0:binw:500e-3;

        % figure,plot(tpsth,psthHS)
        % ylabel('SpikeRate [spikes/s]')
        % xlabel('Time [s]')

        %freq=(0:(numel(tpsth)-1))/(binw)/numel(tpsth);
        %figure,plot(freq,abs(fft(psth))),xlim([0 2500])    

        %% calculate the VS (vector strength)
        %take the spike times vector and then bin depending on the frequency you
        %are testing (i.e., 10ms corresponds to 100Hz) 
        %then calculate the phase for each of these bins (phi=f*t)
        %then calculate a vector with the phase corresponding to the delay for each spike time
        %this becomes a pie with vectors that will rotate for each 10 ms (100Hz)
        %the vector strenght at this frequency will then be the absolute value of the mean of all the
        %calculated vectors
        for F=1:10000;
            phi1=spiketimesHS*F; %the phase value (in cycles)
            vects1=exp(2*pi*j*phi1);
            VSHS(F)=abs(mean(vects1));

            phi2=spiketimesLS*F; %the phase value (in cycles)
            vects2=exp(2*pi*j*phi2);
            VSLS(F)=abs(mean(vects2));
        end

        PAN_HS(:,ch)=psthHS;
        VSAN_HS(:,ch)=VSHS;

        PAN_LS(:,ch)=psthLS;
        VSAN_LS(:,ch)=VSLS;
    end

    eval(['save(''PAN_HS',num2str(WhichF),'_',num2str(L),'.mat'',''PAN_HS'')'])
    eval(['save(''PAN_LS',num2str(WhichF),'_',num2str(L),'.mat'',''PAN_LS'')'])
    eval(['save(''VSAN_HS',num2str(WhichF),'_',num2str(L),'.mat'',''VSAN_HS'')'])
    eval(['save(''VSAN_LS',num2str(WhichF),'_',num2str(L),'.mat'',''VSAN_LS'')'])
end
%toc

