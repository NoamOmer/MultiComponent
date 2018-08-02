function [fVec,fftMat]=plotfft(x, nFFT, fs, winType, unitsStr)
% plotfft(x, nFFT, fs, winType, unitsStr)

N = length(x);
if size(x, 2) == N
    x = x.';
end
nSigs = size(x,2);

if ~exist('nFFT', 'var') || isempty(nFFT)
    nFFT = length(x);
end
if ~exist('fs', 'var') || isempty(fs)
    fs = 1;
    xlabelStr = 'Freq [\pi]';
else
    xlabelStr = 'Freq [Hz]';
end
if ~exist('winType', 'var') || isempty(winType) || strcmp(winType, 'boxcar')
    winMat = ones(N, nSigs);
elseif strcmp(winType, 'hammming')
    winMat = hamming(N);
    winMat = repmat(winMat, 1, nSigs);
elseif strcmp(winType, 'hann')
    winMat = hann(N);
    winMat = repmat(winMat, 1, nSigs);
end
winMat = winMat / min(N, nFFT);
fVec = (-floor(nFFT/2) : ceil(nFFT/2)-1) / nFFT * fs;
fftMat = fftshift( fft(x .* winMat, nFFT, 1), 1);

if ~exist('unitsStr', 'var') || isempty(unitsStr) || strcmpi(unitsStr, 'db')
    absFFT = db( fftMat );
    ylabelStr = 'Power [dB]';
elseif strcmpi(unitsStr, 'linear')
    absFFT = abs( fftMat );
    ylabelStr = 'Magnitude';
end

hold all
for sigNum = 1:nSigs
    plot(fVec, absFFT(:,sigNum), '.-')
end
hold off
grid minor; legend(arrayfun(@(x) sprintf('sig%d', x), 1:nSigs, 'uniformoutput', false))
xlabel(xlabelStr); ylabel(ylabelStr)
fVec=fVec';