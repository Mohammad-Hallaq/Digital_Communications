clear all; clc;
SNR = [-4 0 2:2:8 9:12]; % signal-to-noise ration (dB)
F = 5; % number of frames
M = length(SNR);
Len = 1e3; % number of data bits per frame
ber = zeros(1, M);
for ii = 1 : M
    n = 0; % number of errors
    for jj = 1 : F
        %fprintf('SNR = %f, frame = %d / %d\n', SNR(ii), jj, F);
        b = randsrc(1,Len,[0,1]); % information bits
        t = 2 * b - 1; % BPSK modulation
        r = awgn(t,SNR(ii)); % adding noise
        bb = r > 0; % demodulation and decision
        n = n + nnz(bb-b);
    end
    ber(ii) = n / Len / F;
end
figure(1); semilogy(SNR, ber, 'color', rand(1, 3)); grid on; hold on;
xlabel('SNR (dB)'); ylabel('BER');

%%
clear all; clc;
Er = 200; % number of error needed for each SNR value
maxF = 8e3; % maximum number of transmitted frames for each SNR value
SNR = [-4 0 2:2:8 9:12]; % signal-to-noise ration (dB)
F = 5; % number of frames
M = length(SNR);
Len = 1e3; % number of data bits per frame
ber = zeros(1, M);
for ii = 1 : M
n = 0; % number of errors
for jj = 1 : maxF
fprintf('SNR = %f, frame = %d / %d\n', SNR(ii), jj, maxF);

if n < Er || jj <= maxF
    b = randsrc(1,Len); % information bits
    t = 2 * b - 1; % BPSK modulation
    r = awgn(t,SNR(ii)); % adding noise
    bb = r > 0; % demodulation and decision
    n = n + nnz(bb-b);
end

end
ber(ii) = n / Len / F;
end
figure(1); semilogy(SNR, ber, 'color', rand(1, 3)); grid on; hold on;
xlabel('SNR (dB)'); ylabel('BER');
%%
clear all
clc
SNR = [-4 0 2:2:8 9:12] -3; % signal-to-noise ration (dB)
