close all
clear all

%new version that takes into account the "logical AND" at
%the level of the CN i.e. coincidence detection
%and then puts in the CN and IC model

tic
load Fcnew.mat


FS=100000; %for now with the reshaping
nrep=1;
%ratios for LS/HS in pct
TF=20; %total no of fibers synapsing onto CN
%% should all have odd exponents
rLS=round(15/100*TF); %in percentage of total population
rMS=round(20/100*TF)-1; %-1 done here to keep it odd in case of multiplication
rHS=round(65/100*TF);
TFtot=19; %
sumon=1; %if sumon, then use addition of fibers, else multiplication. If multiplication, only odd numbers of percentages work

gL=[0:10:100];

for r=1:numel(gL);
    LC=gL(r);
    LC
               
    eval(['load ANHS',num2str(LC),'.mat'])
    eval(['load ANMS',num2str(LC),'.mat'])
    eval(['load ANLS',num2str(LC),'.mat'])
    
    %try here to subtract the SR, so that baseline is at the max SR
    ANHS=ANHS-100;
    ANLS=ANLS-100;
    ANMS=ANMS-100;
    
    %300 for the 1kHz resp/4kHz resp 165
    LHS(:,r)=ANHS(:,165); %1kHz responses
    LMS(:,r)=ANMS(:,165); %1kHz responses
    LLS(:,r)=ANLS(:,165); 
    
    Acn=1.5;
    Aic=1;
    Scn=0.6;
    Sic=1.5;
    Dcn=1e-3;
    Dic=2e-3;
    Tex=0.5e-3;
    Tin=2e-3;
    t=[0:size(ANHS)-1]/FS;

    Inhcn=[zeros(1,round(Dcn*FS))  Scn*(1/Tin^2)*(t).*exp(-(t)/Tin)];
    Inhcn(end-round(Dcn*FS)+1:end)=[];

    Inhic=[zeros(1,round(Dic*FS))  Sic*(1/Tin^2)*(t).*exp(-(t)/Tin)];
    Inhic(end-round(Dic*FS)+1:end)=[];

    for n=1:467 %the other sections have not been calculated
        display(num2str(n))
        
        %do first the coincidence detection to sum up LS/HS fibers
        %AN(:,n)=((ANLS(:,n).^rLS).*(ANMS(:,n).^rMS)...
        %    .*(ANHS(:,n).^rHS))./sqrt(rLS+rMS+rHS); 
        
        if sumon==1
             AN=(rLS*ANLS+rHS*ANHS+rMS*ANMS)/(TFtot);
        else
             Part=(ANLS.^rLS).*(ANMS.^rMS).*(ANHS.^rHS);
             AN=sign(Part).*(abs(Part)).^(1/(rLS+rMS+rHS));  
        end
        
        Rcn(:,n)=Acn*(conv((1/Tex^2)*t.*exp(-t/Tex),AN(:,n))-conv(Inhcn,circshift(AN(:,n),round(Dcn*FS))));
        Ric(:,n)=Aic*(conv((1/Tex^2)*t.*exp(-t/Tex),Rcn(:,n))-conv(Inhic,circshift(Rcn(:,n),round(Dic*FS))));
     end

CN4k(:,r)=Rcn(:,165); %1kHz responses
IC4k(:,r)=Ric(:,165);
if r==7
    AN60=AN;
    Rcn60=Rcn;
    Ric60=Ric;
end
    
%% here only summed until CF corresponding to 175Hz
IC(:,r)=sum(Ric(1:size(AN,1),1:430),2); 
W1(:,r)=sum(AN(1:size(AN,1),1:430),2); 
end

save('IC_SMix.mat','Rcn60','Ric60','IC','AN60','W1')
%eval(['save(''responses1k',num2str(LC),'.mat'',''LHS'',''LMS'',''LLS'',''CNLS'',''CNHS'',''CILS'',''CIHS'',''CIsum'')'])
%eval(['save(''IC',num2str(LC),'nn.mat'',''IC'')'])
%eval(['save(''ICLS',num2str(LC),'nn.mat'',''ICLS'')'])
%eval(['save(''ICHS',num2str(LC),'nn.mat'',''ICHS'')'])

%end 