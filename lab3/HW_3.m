clc;
close all;
clear;

M = 1000;
fs = 16;
b = randsrc(1,M,[0, 1]);
t = 2*b -1;
t_ovsample = kron(t, [1 zeros(1,fs-1)]);%upsample
R = [0, 0.5, 1];
alpha = ["alpha = 0", "alpha = 0.5", "alpha = 1"];
N_T = 3;
RATE = fs;
T = 1;
for i=1:3
    RC_fir = rcosfir(R(i),N_T,RATE,T,'normal');
    y = conv(t_ovsample,RC_fir);
    subplot(3,1,i);
    periodogram(y);
    title(alpha(i));
end

T_DELAY = 2;
TOL = 0.03;

for j=1:3
    [NUM, DEN] = rcosiir(R(j), T_DELAY, RATE, T, TOL,'normal');
    y_2 = filter(NUM, DEN, t_ovsample);
    subplot(3,1,j);
    periodogram(y_2);
    title(alpha(j));
    
    
end
%%
M = 1000;
fs = 16;
b = randsrc(1,M,[0, 1]);
t = 2*b -1;
t_ovsample = kron(t, [1 zeros(1,fs-1)]);
R = 1;
N_T = 3;
RATE = fs;
T = 1;
RC_fir = rcosfir(R,N_T,RATE,T,'normal');
plot(RC_fir);
seq = conv(t_ovsample,RC_fir);
t =0:15;
for k=4:length(b)-4
    
     
    plot(t, seq(9+16*k:24+k*16), 'black')
    hold on
   
       
    
end
%%
clc;
T=1; 
RATE=16; 
N_T=6;  
R=[0 0.5 1]; 
LenFrame=100000; 
ebno=0:1:8;
Sampling=[1 zeros(1,RATE-1)]; 
for k=1:length(R)
    Flter=rcosfir(R(k),N_T,RATE,T,'sqrt');
    for j=1:length(ebno) 
            NumE=0;
            NumFramw=0;
            while (NumFramw<1000)&&(NumE<200)   
                 optemuminstance=1+RATE*(N_T)*2; 
                 Signal=randsrc(1,LenFrame,[-1 1]);  
                 SampledSignal=kron(Signal,Sampling);   
                 ConvSignal=conv(Flter,SampledSignal);  
                 SignalNoised=awgn(ConvSignal,ebno(j) + 3); 
                 RecivedSig=conv(Flter,SignalNoised);        
                 receive_signal=RecivedSig(optemuminstance:RATE:length(RecivedSig)-optemuminstance);
                 receive_element=(receive_signal>0)-(receive_signal<0);
                 NumE=NumE+nnz(receive_element-Signal);
                 NumFramw=NumFramw+1;
            end
            ber(k,j)=NumE/(NumFramw*LenFrame);
    end
end

for j=1:length(R)
    semilogy(ebno,ber);
    hold on
end
y=berawgn(ebno ,'Pam',2);
semilogy(ebno,y)
title('BIT Error Rate');
xlabel('EB/N0');
ylabel('Pe');
legend('theoritical','beta=1','beta=0.5','beta=0')
%%
T=1; 
A=1; 
sps=16; 
n_t=6; 
beta=[0 0.5 1]; 
LenFrame=100000; 
ebno=0:1:8;  
delay=2;
Sampling=[1 zeros(1,sps-1)];
for k=1:length(beta)
    Flter=rcosfir(beta(k),n_t,sps,T,'sqrt');
    for j=1:length(ebno) 
            NumE=0;
            NumFramw=0;
            while (NumFramw<1000)&&(NumE<200)   
                 optemuminstance=delay+1+sps*(n_t)*2; 
                 Signal=randsrc(1,LenFrame,[-1 1]);  
                 SampledSignal=kron(Signal,Sampling);   
                 ConvSignal=conv(Flter,SampledSignal);  
                 SignalNoised=awgn(ConvSignal,ebno(j) + 3); 
                 RecivedSig=conv(Flter,SignalNoised);        
                 receive_signal=RecivedSig(optemuminstance:sps:length(RecivedSig)-optemuminstance);
                 receive_element=(receive_signal>0)-(receive_signal<0);
                 NumE=NumE+nnz(receive_element-Signal);
                 NumFramw=NumFramw+1;
            end
            ber(k,j)=NumE/(NumFramw*LenFrame);
    end
end

for j=1:length(beta)
    semilogy(ebno,ber);
    hold on
end
y=berawgn(ebno ,'Pam',2);
semilogy(ebno,y)
title('BIT Error Rate');
xlabel('EB/N0');
ylabel('Pe');
legend('theoritical','beta=1','beta=0.5','beta=0')
