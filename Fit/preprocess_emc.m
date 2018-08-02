function [E,E_noisy]=preprocess_emc(axisT2,emc,firstTE,snr,denoise_flag,ETL,type)

E=[]; E_noisy=[];

% Normalize EMC databases
% ----------------------------------------------------------------------------------------
emcNorm = normalize_5D_EMC(emc);


% Calculate a decay matrix for each simulated EMC vector accorting to the first TE (instead of the EMC starting at '1')
% ----------------------------------------------------------------------------------------
decayVec        = transpose(exp(-firstTE./(1e3*axisT2)));
decayMat        = repmat(decayVec,1,ETL);
emcNormDecay = squeeze(emcNorm).*decayMat;
%emcNormDecay = emcNorm ; % a test remember to delete this row ! 

% Generate output EMC operator (Aeq)
% ----------------------------------------------------------------------------------------
if strcmp(type,'exprmt')
    % Create a noisy EMC DB (and de-noise) for constructing an "experimental" emc
    [emcNormDecayNoisy, emcNormDecayDeNoised] = noise_and_denoise_EMC_DB(snr,emcNormDecay,denoise_flag,ETL);
    E_noisy = (squeeze(emcNormDecayNoisy)); % Verify with Noam
	E		= emcNormDecay;
elseif strcmp(type,'forFit')
    E       = emcNormDecay;
end


end