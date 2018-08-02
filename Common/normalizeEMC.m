% Normalize EMC(B0,B1,T1,T2,TE)

% The following operation factors out both the baseline amplitude and baseline phase
% leaving the first elements equal to 1.0, while the subsequent ones remain complex

function emcNorm = normalizeEMC(emc)

emcNorm=emc./repmat(emc(:,1),1,size(emc,2));

end

