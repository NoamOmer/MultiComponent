% Normalize EMC(B0,B1,T1,T2,TE)

% The following operation factors out both the baseline amplitude and baseline phase
% leaving the first elements equal to 1.0, while the subsequent ones remain complex

function sim_mat = normalize_4D_EMC(sim_mat)
for B0_idx = 1:size(sim_mat,1)
for B1_idx = 1:size(sim_mat,2)
for T2_idx = 1:size(sim_mat,3)
	tmp = squeeze(sim_mat(B0_idx,B1_idx,T2_idx,:));
	tmp = tmp / tmp(1);
	if (sum(tmp<0))
		error('EMC error');
	end;
	sim_mat(B0_idx,B1_idx,T2_idx,:) = tmp;
end;
end;
end;