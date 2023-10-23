clear all; clc;
M = 1e3 ;
Fs = 10;
T =1;
ber = zeros(1, M);
b = randsrc(1,M,[0,1]);
t = 2 * b - 1;
lc = kron(t, ones(1,Fs));

spec_lc = fftshift(fft(lc))/Fs;

ESD = abs(spec_lc).*abs(spec_lc);
PSD = ESD/T;

%plot(ESD);

K = 100;
%PSD_vec = zeros(K,1e4);
PSD_vec = zeros(1,1e4);
for i = 1 : K 
    
    b = randsrc(1,M,[0,1]);
    t = 2 * b - 1;
    lc = kron(t, ones(1,Fs));
    spec_lc = fftshift(fft(lc))/Fs;
    
    ESD = abs(spec_lc).*abs(spec_lc);
    
    PSD_vec = PSD_vec + ESD/1000/T;
    
    
end

PSD_est = PSD_vec/K;

sym = 1:length(PSD_vec);

f = (sym-length(PSD_vec)/2)/(T*1000);

PSD_theo = 1*T*(sinc(f*T)).^2;

figure;
hold on; grid on;
plot(f,PSD_est);
plot(f, PSD_theo);
xlabel("Frequency"); ylabel("PSD"); title("PSD For NRZ-L");
legend("Estimated", "Theoretical");

%%
clear all; clc;
M = 1e3 ;
Fs = 10;
T =1;
ber = zeros(1, M);
b = randsrc(1,M,[0,1]);
t = 2 * b - 1;
lc = kron(t, [1,1,1,1,1,0,0,0,0,0]);
spec_lc = fftshift(fft(lc))/Fs;

ESD = abs(spec_lc).*abs(spec_lc);
PSD = ESD/T;

%plot(ESD);

K = 100;
PSD_vec = zeros(K,1e4);

for i = 1 : K 
    
    b = randsrc(1,M,[0,1]);
    t = 2 * b - 1;
    lc = kron(t, [1,1,1,1,1,0,0,0,0,0]);
    spec_lc = fftshift(fft(lc))/Fs;
    
    ESD = abs(spec_lc).*abs(spec_lc);
    
    PSD_vec(i,:) = ESD/1000/T;
    
    
end


PSD_est = mean(PSD_vec);

sym = 1:length(PSD_vec);

f = (sym-length(PSD_vec)/2)/(T*1000);
PSD_theo = (1/4)*T*(sinc(f*T/2)).^2;

figure;
plot(f, PSD_est);
hold on; grid on;
plot(f, PSD_theo);
xlabel("Frequency"); ylabel("PSD"); title("PSD For Polar RZ");
legend("Estimated", "Theoretical");

%%
clear all; clc;
M = 1e3 ;
Fs = 10;
T =1;
ber = zeros(1, M);
b = randsrc(1,M,[0,1]);
t = 2 * b - 1;
lc = kron(t, [1,1,1,1,1,-1,-1,-1,-1,-1]);

spec_lc = fftshift(fft(lc))/Fs;

ESD = abs(spec_lc).*abs(spec_lc);
PSD = ESD/T;

K = 100;
PSD_vec = zeros(K,1e4);

for i = 1 : K 
    
    b = randsrc(1,M,[0,1]);
    t = 2 * b - 1;
    lc = kron(t, [1,1,1,1,1,-1,-1,-1,-1,-1]);
    spec_lc = fftshift(fft(lc))/Fs;
    
    ESD = abs(spec_lc).*abs(spec_lc);
    
    PSD_vec(i,:) = ESD/1000/T;
    
    
end


PSD_est = mean(PSD_vec);

sym = 1:length(PSD_vec);

f = (sym-length(PSD_vec)/2)/(T*1000);
PSD_theo = 1*T*(sinc(f*T/2).^2).* (sin(pi*f*T/2)).^2;

figure;
hold on;
plot(f, PSD_est);
hold on; grid on;
plot(f, PSD_theo);
xlabel("Frequency"); ylabel("PSD"); title("PSD For Manchester");
legend("Estimated", "Theoretical");

%%
clc
close all

EbNo = 0:12;
X = 10.^(EbNo/10);
BER = qfunc(sqrt(2*X));

figure(1); semilogy(EbNo, BER, 'color', rand(1, 3)); grid on; hold on;
xlabel('Eb/N0(dB)'); ylabel('BER');title("NRZ BER");
%%
clc
close all

EbNo = 0:12;
X = 10.^(EbNo/10);
BER = qfunc(sqrt(2*X));

figure(1); semilogy(EbNo, BER, 'color', rand(1, 3)); grid on; hold on;
xlabel('Eb/N0(dB)'); ylabel('BER');title("POLAR-RZ BER");
%%
clc
close all

EbNo = 0:12;
X = 10.^(EbNo/10);
BER = qfunc(sqrt(2*X));

figure(1); semilogy(EbNo, BER, 'color', rand(1, 3)); grid on; hold on;
xlabel('Eb/N0(dB)'); ylabel('BER');title("MANCHESTER BER");
