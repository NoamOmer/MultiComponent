function params=getParams()

	simData = '2016_09_22_T2_1_1000_B1_100_for_MultiT2Fit_400um_L_20mm_ABS';
	
	% Synthesized (experimental) EMC params
    % ----------------------------------------------------------------------------------------
    params.nExprT2                  = 1;					% # of T2 components in synthesized vector
    params.peakWidth                = 20;		
    params.mixMode                  = 2;					% [1: Series of single T2 values | 2: Series of gaussians distribution around each T2]
    params.loc					    = [ 600   ,0  ,0	];		% [ms]		T2 components at location 'loc'
    params.relativeFract		    = [ 1	 ,0  ,0		];	% [none]	relative fraction of each component
    params.winWidth					= params.peakWidth*ones(1,params.nExprT2);	 % [none]	width of distribution around each T2 component
    params.neighbor_mrgn            = 0;					% [%] variability of T2 value and fraction between neighboring pixels (20)
    params.n_neighbors              = 0;					% Number of neighboring pixels (0: no neighbors) - optional 8
    params.snr                      = 0;					% 0: no noise. MRI typical SNR 30-100.
    params.denoise_flag.denoise_op  = 0;					% [0 | 1 | 2 : which denoising algorithm to use]
    params.denoise_flag.r           = 5;
    params.isExperimental           = 1;
    % ----------------------------------------------------------------------------------------
    
    

    % Fitting operator params
    % ----------------------------------------------------------------------------------------
    params.d_T2                     = 1;
    params.maxT2                    = 1000;
    params.minT2                    = 1;	
    params.sparsOrder               = 3;
    params.min_T2_distance          = 2;
    % ----------------------------------------------------------------------------------------
    

	% Solver params
    % ----------------------------------------------------------------------------------------
    
	% Assuming averaging, SNR should be x3 for 2D data, or x3*sqrt(3) for 3D data
    % Don't use this though -- it's unclear how to denoise a noisy experimental EMC.
  
	params.Sum_fitted_w_Eql_to_one_flag   = 1;     % [0 | 1] Fitting constraint: set sum of fitted weights == 1
	params.Tikh_reg_flag                  = 0;     % [0 | 1] Enable / disable the use of Tikhonov regularization in fitting
	params.lambda_Tikh                    = 1e-6;  % Tikhonov regularization weight
	params.w0_Tikh                        = 0;     % Apply Tikhonov around non zero baseline value
	params.lambda_L1                      = 1e-5;  % L1       regularization weight
	params.nT2s_to_fit_to                 = 4;     % [0: all T2s in the EMC DB  |  N: N*ETL T2s, spread logarithmically]
	params.fit_TYPE                       = 2;     % [1: Linear  ||  2: Quadratic  || 3: OMP]
	params.fit_nRep                       = 1;     % Repeat fit in case it can use prior info from previous runs 
												   % quadprog with nRep=2 seems to de-noise the vector of weights
	% ----------------------------------------------------------------------------------------
end




%
% 			nExprT2					 = length(loc);
%             switch nExprT2
%                 case 1
%                     % loc          = round(([0.5 ]*maxT2 - minT2+1)*(n_sprs_T2s/maxT2));
%                     loc            = [60];
%                     relativeFract  = [1];
%                     winWidth       = [40];
%
%                 case 2
%                     loc            = [30 , 64];
%                     relativeFract  = [0.5, 0.5];
%                     winWidth       = round([peakWidth, peakWidth]); %*(n_sprs_T2s/maxT2));
%
%                 case 3
%                     loc            = [30   ,64 , 800 ];   % [ms]
%                     relativeFract  = [0.10 ,0.75, 0.15];
%                     winWidth       = round([peakWidth*0.5 ,peakWidth*2,peakWidth*12]);
%
%                 case 4
%                     loc            = [15   ,80 ,200, 500];   % [ms]
%                     relativeFract  = [0.6 ,0.1 ,0.2 ,0.1];
%                     winWidth       = [30  ,20  ,30  ,20 ];
% 			end
