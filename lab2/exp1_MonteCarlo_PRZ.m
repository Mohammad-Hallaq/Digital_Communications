function [ber, numBits] = exp1_MonteCarlo(EbNo, maxNumErrs, maxNumBits)
% Import Java class for BERTool.
import com.mathworks.toolbox.comm.BERTool;

% Initialize variables related to exit criteria.
berVec = zeros(3,1); % Updated BER values
M = 1e3;
Fs = 10;
g = [sqrt(2)*ones(1,Fs/2) zeros(1, Fs/2)]/sqrt(Fs);
% --- Set up parameters. ---
% --- INSERT YOUR CODE HERE.
% Simulate until number of errors exceeds maxNumErrs
% or number of bits processed exceeds maxNumBits.
while((berVec(2) < maxNumErrs) && (berVec(3) < maxNumBits))
        b = randsrc(1,M,[0,1]);
        t = 2 * b - 1;
        lc = sqrt(2)*kron(t, [1,1,1,1,1,0,0,0,0,0])/sqrt(Fs);
        r = awgn(lc,EbNo + 3);
        rr = reshape(r,Fs,M);
        res = g*rr;
        bb = res > 0; % demodulation and decision
        berVec(2) = berVec(2) + nnz(bb-b);
        berVec(3) = berVec(3) + M;
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