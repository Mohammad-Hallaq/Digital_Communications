function [ber, numBits] = exp1_MonteCarlo(EbNo, maxNumErrs, maxNumBits)
% Import Java class for BERTool.
import com.mathworks.toolbox.comm.BERTool;

% Initialize variables related to exit criteria.
berVec = zeros(3,1); % Updated BER values
Len = 1e3;
fs = 16;
R = 0.3;
N_T = 3;
RATE = fs;
T = 1;
g = ones(1,fs)/sqrt(fs);
% --- Set up parameters. ---
% --- INSERT YOUR CODE HERE.
% Simulate until number of errors exceeds maxNumErrs
% or number of bits processed exceeds maxNumBits.
while((berVec(2) < maxNumErrs) && (berVec(3) < maxNumBits))
        b = randsrc(1,Len,[0,1]); % information bits
        t = 2 * b - 1; % BPSK modulation
        t_ovsample = upsample(t, fs);
        RC_fir = rcosfir(R,N_T,RATE,T,'sqrt');
        seq = conv(t_ovsample,RC_fir);
        r = awgn(seq,EbNo + 3); % adding noise
        RC_fir_2 = rcosfir(R,N_T,RATE,T,'sqrt');
        out =  conv(r,RC_fir_2);
        rr = reshape(out,fs,[]);
        res = g*rr;
        bb = res > 0;
        berVec(2) = berVec(2) + nnz(bb-b);
        berVec(3) = berVec(3) + Len;
   % Check if the user clicked the Stop button of BERTool.
   if (BERTool.getSimulationStop)
      break;
   end

   % --- Proceed with simulation.
   % --- Be sure to update totErr and numBits.
   % --- INSERT YOUR CODE HERE.
end % End of loop

berVec(1) = berVec(2)/berVec(3);

% Assign values to the output variables.
ber = berVec(1);
numBits = berVec(3);