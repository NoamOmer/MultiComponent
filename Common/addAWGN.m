function noisySig = addAWGN(signal,desierdSNR)
% ---------------------------------------------------------------------
% Function Description :
% This function apply white Gaussian noise given a specified SNR value 
% assuming input signal is 1D and real valued
%
% Important Note :
% SD    = sqrt(1/n*sum((X-XMean)^2));
% RMS   = sqrt(1/n*sum(X^2)); 
% For   XMean = 0 (which is the case for white Gaussian noise) -> SD = RMS 
%
% Writen by Noam Omer 
% ---------------------------------------------------------------------

sigLen      = length(signal);                   % Length
noiseVec    = randn(size(signal));              % Orignal noise
pwrSig      = sqrt(sum(signal.^2))/sigLen;      % Signal power
pwrNoise    = sqrt(sum(noiseVec.^2))/sigLen;    % Noise power (RMS=STD)

if desierdSNR ~= 0
    scaleFactor = (pwrSig/pwrNoise)/desierdSNR; % Find scale factor
    noiseVec = scaleFactor*noiseVec;            % Scale the noise
    noisySig = signal + noiseVec;               % Add noise
else
    noisySig = noiseVec;                       % Noise only
end

return;