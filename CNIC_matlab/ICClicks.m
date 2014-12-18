function ICClicks(name,Neuro)

%single channel load.
L=0:10:100;
tic
load CFs.mat
%L=75; 
%name='NHClicksME';
%Neuro=1;
%Neuro=3;
Ncut=find(CF>1000,1,'last') %indices below this boundary have HF neuropathy

TF=20; %total no of fibers synapsing onto CN
if Neuro==1
    rLS=0;%round(15/100*TF); %in percentage of total population
    rMS=0;%round(20/100*TF)-1;
    rHS=round(65/100*TF);
elseif Neuro==2
    rLS=0;%round(15/100*TF); %in percentage of total population
    rMS=0;%round(20/100*TF)-1;
    rHS=round(round(65/100*TF)/2);
elseif Neuro==3 || 4  %freq specific neuropathy
    rLSN=0;
    rMSN=0;
    rHSN=0;
    rLS=round(15/100*TF); %in percentage of total population
    rMS=round(20/100*TF)-1;
    rHS=round(65/100*TF);    
else
    rLS=round(15/100*TF); %in percentage of total population
    rMS=round(20/100*TF)-1;
    rHS=round(65/100*TF);    
end
TFtot=19;
%ratios for LS/HS in pct

IC=0; W1=0; CN=0;


for k=1:numel(L)
    disp(num2str(L(k)))
     
    load (['../ANerve_matlab/out/Clicks/',name,'ANLS_',num2str(L(k)),'.mat'])
    load (['../ANerve_matlab/out/Clicks/',name,'ANMS_',num2str(L(k)),'.mat'])
    load (['../ANerve_matlab/out/Clicks/',name,'ANHS_',num2str(L(k)),'.mat'])

    
    for n=1:size(HS,2)
        disp(num2str(n))
        %sampling rate down to 20kHz
        ANHS=resample(HS(:,n),1,5);
        ANMS=resample(MS(:,n),1,5);
        ANLS=resample(LS(:,n),1,5);
    
        FS=20000;
       
        Acn=1.5;
        Aic=1;
        Scn=0.6;
        Sic=1.5;
        Dcn=1e-3;
        Dic=2e-3;
        Tex=0.5e-3;
        Tin=2e-3;
        t=[0:size(ANHS,1)-1]/FS';

        Inhcn=[zeros(1,round(Dcn*FS))  Scn*(1/Tin^2)*(t).*exp(-(t)/Tin)];
        Inhcn(end-round(Dcn*FS)+1:end)=[];

        Inhic=[zeros(1,round(Dic*FS))  Sic*(1/Tin^2)*(t).*exp(-(t)/Tin)];
        Inhic(end-round(Dic*FS)+1:end)=[];
         
        if Neuro==3
            if  CF(n)>Ncut %for the neuropathy region HF
                AN=(rLSN*ANLS+rHSN*ANHS+rMSN*ANMS)/(TFtot);
            else
                AN=(rLS*ANLS+rHS*ANHS+rMS*ANMS)/(TFtot);
            end
        elseif Neuro==4
           	if  CF(n)<Ncut %for the neuropathy region LF
                AN=(rLSN*ANLS+rHSN*ANHS+rMSN*ANMS)/(TFtot);
            else
                AN=(rLS*ANLS+rHS*ANHS+rMS*ANMS)/(TFtot);
            end     
        else
            AN=(rLS*ANLS+rHS*ANHS+rMS*ANMS)/(TFtot);
        end
        
        Rcn=Acn*(conv((1/Tex^2)*t.*exp(-t/Tex),AN)-conv(Inhcn,circshift(AN,round(Dcn*FS))));
        Ric=Aic*(conv((1/Tex^2)*t.*exp(-t/Tex),Rcn)-conv(Inhic,circshift(Rcn,round(Dic*FS))));

    %% here only summed until 175Hz
    %population responses
    if n<=433 %the summed response for the case
        W1=W1+AN(1:numel(AN)); %add them up one by one
        CN=CN+Rcn(1:numel(AN)); %add them up one by one
        IC=IC+Ric(1:numel(AN)); %add them up one by one
    end

     RicF(n,:)=Ric(1:numel(AN));
     RcnF(n,:)=Rcn(1:numel(AN));
    
end
toc

if Neuro==1
    save(['./out/Clicks/IC',name,num2str(L(k)),'_HS.mat'],'IC','CN','W1')
    save(['./out/Clicks/IC',name,num2str(L(k)),'resp_HS.mat'],'RicF','RcnF')
elseif Neuro==2
    save(['./out/Clicks/IC',name,num2str(L(k)),'_HS50.mat'],'IC','CN','W1')
    save(['./out/Clicks/IC',name,num2str(L(k)),'resp_HS50.mat'],'RicF','RcnF')
elseif Neuro==3
    save(['./out/Clicks/IC',name,num2str(L(k)),'_HSHF.mat'],'IC','CN','W1')
    save(['./out/Clicks/IC',name,num2str(L(k)),'resp_HSHF.mat'],'RicF','RcnF')
elseif Neuro==4
    save(['./out/Clicks/IC',name,num2str(L(k)),'_HSLF.mat'],'IC','CN','W1')
    save(['./out/Clicks/IC',name,num2str(L(k)),'resp_HSLF.mat'],'RicF','RcnF')
else
    save(['./out/Clicks/IC',name,num2str(L(k)),'_Mix.mat'],'IC','CN','W1')
    save(['./out/Clicks/IC',name,num2str(L(k)),'resp_Mix.mat'],'RicF','RcnF')
end

end %end of all levels