function [emc_db_wNoise,emc_db_DeNoised] = noise_and_denoise_EMC_DB(snr,emc_db,denoise_opts,inETL,plot_f)

if (snr == 0)
	emc_db_wNoise   = emc_db;
	emc_db_DeNoised = emc_db;
	return;
end;

if (ndims(emc_db) == 2)
	tmp(1,1,1,:,:) = emc_db;
	emc_db = tmp;
end;
[nB0,nB1,nT1s,nT2s,ETL] = size(emc_db);
if (inETL ~= ETL) || (nB0 ~= 1)
	error('Dimensions mismatch or unsupported functionality');
end;

global reset_RAND_flag;
if exist('reset_RAND_flag','var') && ~isempty(reset_RAND_flag) && (reset_RAND_flag == 1)
rng('default');  % reset random number generator
reset_RAND_flag = 0;
end;

for B0_idx = 1
for B1_idx = 1:nB1	
	emc_db_T1T2_wNoise = zeros([nT1s,nT2s,ETL]);
	for T1_idx = 1:nT1s
	for T2_idx = 1:nT2s
		emc_vec  = transpose(squeeze(emc_db(B0_idx,B1_idx,T1_idx,T2_idx,:)));
		noiseVec = transpose(wgn(ETL,1,1,'real'));			% Generate a noise vector
		noiseVec = noiseVec/max(abs(noiseVec));				% Normalize the noise vector to 0...1
		noiseVec = noiseVec * (emc_vec(1)/2);				% Scale the noise vector to the current EMC
		emc_db_T1T2_wNoise(T1_idx,T2_idx,:) = emc_vec + noiseVec / snr;
	end;
	end;
% 	figure(1); imagesc(squeeze(emc_db_T1T2_wNoise(1,:,:))- squeeze(emc_db(1,1,1,:,:))); axis square; colorbar; %pause(1);

	% De-noise
	switch denoise_opts.denoise_op
	case 0
		emc_db_T1T2_DeNoised = emc_db_T1T2_wNoise;

	case 1
		error('This code probably requires 2D matrices and thus another loop around T1 index');
		Xtil = emc_db_T1T2_wNoise;
		% r: Estimate of (effective) rank of signal matrix
		r    = denoise_opts.r;								% r=5 produces best results
		[Shat,relmse_hat,mse_hat] = optshrink(Xtil,r);
		emc_db_T1T2_DeNoised = Shat;

	case 2
		error('This code probably requires 2D matrices and thus another loop around T1 index');
		tmp_mat(1,:,:) = emc_db_T1T2_wNoise;
		param.E      = 1;
		param.d      = tmp_mat;
		param.lambda = 1e-2;         % [0.001 ... 0.01] ?
		param.nite   = 40;
		param.stopping_criteria_threshold = 1e-10;

		tmp_mat = lr_mat_ist(param);
		emc_db_T1T2_DeNoised = squeeze(tmp_mat);
	end;

	% Set back into output parameters
	emc_db_wNoise  (B0_idx,B1_idx,:,:,:) = emc_db_T1T2_wNoise;
	emc_db_DeNoised(B0_idx,B1_idx,:,:,:) = emc_db_T1T2_DeNoised;
	
	if (exist('plot_f','var') && plot_f)
		emc_db_T1T2 = emc_db(B0_idx,B1_idx,:,:,:);
		T1idx=ceil(nT1s/2);
		figure;
		subplot(221); imagesc(squeeze(emc_db_T1T2         (T1idx,:,:))); colorbar; axis square; title('EMC DB');
		subplot(222); imagesc(squeeze(emc_db_T1T2_wNoise  (T1idx,:,:))); colorbar; axis square; title('EMC DB w Noise');
		subplot(223); imagesc(squeeze(emc_db_T1T2_DeNoised(T1idx,:,:))); colorbar; axis square; title('EMC DB de-Noised');
		subplot(224); imagesc(squeeze(emc_db_T1T2_DeNoised(T1idx,:,:)) - ...
							  squeeze(emc_db_T1T2         (T1idx,:,:))); colorbar; axis square; title('de-Noised - Orig');
	end;

end; % for B1_idx
end; % for B0_idx

return;




