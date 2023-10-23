clc;
clear;
close all;
M = 1000;
osf = 16;
T = 1;
W = 1/(2*T);
b = randsrc(1,M,[0, 1]);
s = 2*b -1;
s_ovsample = upsample(s, osf);
t = -5*T:1/osf:(5*T - 1/osf);
x = sinc(2*W*t) + sinc(2*W*t - T);
figure;
plot(t, x)
grid on;

filtered_s = conv(s_ovsample,x);

range = 0:osf-1;
for k=3:length(b)-3
    
    

     
    plot(range, filtered_s(9+16*k:8+16+k*16))
    hold on
   
       
    
end
%%
clc;
clear;
close all;
M = 1000;
osf = 16;
T = 1;
W = 1/(2*T);
b = randsrc(1,M,[0, 1]);

pr_co = zeros(1, length(b));
for i=1:length(b)
    if i == 1 
    pr_co(i) = xor(0, b(i));
   
    else
    pr_co(i) = xor(pr_co(i-1), b(i));
    end
end

a = 2*pr_co - 1;

t = -5*T:1/osf:(5*T - 1/osf);
x = sinc(2*W*t) + sinc(2*W*t - T);

x_2 = ifft((abs(fft(x))).^0.5);
figure;
plot(t, x_2)
grid on;

a_ovsample = upsample(a, osf);
filtered_a = conv(a_ovsample,x_2);


