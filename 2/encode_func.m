function [ series ] = encode_func( input_vec)

[H,G,n,k] = hammgen(3);

t=length(input_vec)/k;

series=0;
for i=1:t
 add1 = rem(input_vec(1+(i-1)*k:i*k) * G,2);
 series=[series add1];
end
series=series(2:length(series));
end

