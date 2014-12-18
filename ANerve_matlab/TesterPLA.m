close all
clear all
FS=10000;
cf=1000;

A=linspace(4000,40,3000);
B=linspace(40,60,300);
%synout=[60*ones(100,1) ; 4000*ones(2000,1) ; A'.*ones(3000,1) ;  B'.*ones(300,1) ; 60*ones(10000,1)];
synout=[0 ; ones(1000,1) ; zeros(100,1)];

sampFreq=100000;
binwidth = 1/sampFreq; %should be new Fs if it is available

%synout=resample(synout,sampFreq,FS);

%alpha1 = 5e-6*100e3; beta1 = 5e-4; I1 = 0;*/ /* older version, 2012 and before */
alpha1 = 2.5e-6*FS; beta1 = 5e-4; I1 = 0; %old rate
alpha2 = 1e-2*FS; beta2 = 1e-1; I2 = 0;        
delaypoint=floor(7500/(cf/1e3));

I1=0; %at time zero the integrator is zero
I2=0;

for n=1:numel(synout)
    disp(num2str(n))
    sout1(n)=synout(n)-alpha1*I1;
        if sout1(n)<0
            sout1(n)=0;
        end
    sout2(n)=synout(n)-alpha2*I2;
        if sout2(n)<0
            sout2(n)=0;
        end
    I1=0; I2=0;
    for k=1:n
       I1=I1+sout1(k)*binwidth/((n-k)*binwidth+beta1); %Drew and Abbott PLA function
       I2=I2+sout2(k)*binwidth/((n-k)*binwidth+beta2); 
    end
        sout(n)=sout1(n)+sout2(n); %((sout1(n)+sout2(n)));
end

figure,plot(synout,'linew',2),hold on,plot(sout,'r--'),plot(sout1,'k--'),plot(sout2,'m:'),plot(sout/2,'r','linew',2)
legend('input','Zilanymethod=PLA1+PLA2','PLA1','PLA2','better:(PLA1+PLA2)/2')

