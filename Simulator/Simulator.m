classdef Simulator
    %% Description:
    %       This fucntion generates an artificial EMC, "acquired" for a sample containing a distribution of T2 values
    %       (single values, gaussian distribution etc.). This emc is added with noise according to input SNR.
    %       Specific calculation of this emc is done by geneerating a noisy EMC DB, and multiplying it by a chosen
    %       set of T2 components
    % ----------------------------------------------------------------------------------------
    % ----------------------------------------------------------------------------------------
    
    %% Static functions required:
	%		setParams()
    %       sim_mat = normalize_5D_EMC(sim_mat)
    %       T2_axis = calc_equispaced_T2(T2l,T2h,T2base,nT2s,min_T2_distance,round_f,plot_f)
    %       [emc_db_wNoise,emc_db_DeNoised] = noise_and_denoise_EMC_DB(snr,emc_db,denoise_opts,inETL,plot_f)
    %% ----------------------------------------------------------------------------------------
    %% ----------------------------------------------------------------------------------------
    
    properties
        
        operator = struct(...
							'axisT2',				[],...
							'n_T2s',                [],...
							'd_T2',                 [],...
                            'minT2',                [],...
                            'maxT2',                [],...
                            'sparsOrder',   		[],...
                            'min_T2_distance',		[],...
							'E',                    [] ...
						  );
              
        experimental = struct(...
							'peakWidth',		[],...
							'loc',				[],...
							'relativeFract',	[],...
							'winWidth',			[],...
							'n_Neighbors',		[],...
							'neighbor_mrgn',	[],...
							'snr',				[],...
							'denoise_flag',		[],...
							'mixMode',			[],...
							'nExprT2',			[],...
                            'E_noisy',			[],...
                            'isExperimental',	[] ...
                           );

        DBprop = struct(...
							'firstTE',			[],...            
							'n_B0',				[],...
							'n_B1',				[],...
							'n_T1',				[],...
							'n_T2',				[],...
							'ETL',				[] ...
						); 
    end
    
    methods    

        function obj=loadSimData(obj)
			
			% ----------------------------------------------------------------------
			% Load EMC data from a saved database (i.e. a '.mat' file)
			% ----------------------------------------------------------------------            			
            simData						= load('C:\Users\me\Google Drive\Noam Omer\MultiCompT2\DB\2016_09_22_T2_1_1000_B1_100_for_MultiT2Fit_400um_L_20mm_ABS.mat');      % Nz=4000, EMC_Abs_flag=1, TE=12, ETL=40, RefAng=120, BW=200          
			axis_T2_full				= 1e3*simData.T2_tse_arr;	% converison from mSec to Sec	
            emcT2						= abs(simData.echo_train_modulation);
            [n_B0,n_B1,n_T1,n_T2,ETL]	= size(emcT2);
			firstTE						= 1e3*simData.TE_arr(1);	% converison from mSec to Sec
			B1_1						= find(simData.B1_scaling_arr==1); 
			TE_arr						= simData.TE_arr;
			
			if  isempty(B1_1);  error('invalid B1 scale');    end;
			if  (n_B1 > 1),   emcT2 = emcT2(1,B1_1,:,:,:);  end;
            if  (n_T1 > 1),   emcT2 = emcT2(1,: ,1,:,:);    end;
			
            
			% Normelize original EMC matrix
			% [nbe] cancel preprocess function
			[emcFullNorm,~]=preprocess_emc(axis_T2_full,emcT2,firstTE,[],[],ETL,'forFit');
			
			%----------------------------------------------------------------------
			% Store data in the current instance : 
			%----------------------------------------------------------------------
			obj.DBprop.simData		= simData;            
			obj.DBprop.firstTE		= firstTE;
			obj.DBprop.emcT2_full	= emcT2;
            obj.DBprop.axis_T2_full = axis_T2_full;
			obj.DBprop.emcFullNorm  = emcFullNorm;
			obj.DBprop.TE_arr		= TE_arr;
			%----------------------------------------------------------------------

        end
		function obj=getParams(obj,type)
			
			% Load a set of parameters
			% ----------------------------------------------------------------------------------------
			params=setParams();
			
			
			% Assign each parameter in relevant structure
			% ----------------------------------------------------------------------------------------
			
			% TRUE -> Params for the experimental EMC  (e):
			if strcmp(type,'exprmt')
				
				obj.experimental.peakWidth				 = params.peakWidth;
				obj.experimental.nExprT2				 = params.nExprT2;
				obj.experimental.n_Neighbors			 = params.n_neighbors;
				obj.experimental.neighbor_mrgne			 = params.neighbor_mrgn;
				obj.experimental.snr					 = params.snr;
				obj.experimental.mixMode				 = params.mixMode;
				obj.experimental.denoise_flag.denoise_op = params.denoise_flag.denoise_op;
				obj.experimental.denoise_flag.r			 = params.denoise_flag.r;
				obj.experimental.loc					 = params.loc;
				obj.experimental.relativeFract			 = params.relativeFract;
				obj.experimental.winWidth				 = params.winWidth;
				obj.experimental.nExprT2                 = params.loc;
				
			% FALSE -> Params for the fitting operator (E):
			elseif strcmp(type,'forFit')
				obj.operator.d_T2						 = params.d_T2;
				obj.operator.maxT2						 = params.maxT2;
				obj.operator.minT2						 = params.minT2;
				obj.operator.sparsOrder					 = params.sparsOrder;
				obj.operator.min_T2_distance			 = params.min_T2_distance;

            end
		end
        function obj=genrtNoiseVec(obj)
            
			n_neighbors              = obj.experimental.n_Neighbors;
			neighbor_mrgn            = obj.experimental.neighbor_mrgne;

			neighbor_var = [];
			if (n_neighbors > 0);
				
				% rng('default');                                               % Reset random number generator
				neighbor_var        = transpose(wgn(n_neighbors,1,1,'real'));   % Generate noise vector for 8 neighbors
				neighbor_var        = neighbor_var/(max(abs(neighbor_var)));    % Scale values to -1...+1
				neighbor_var        = neighbor_var*neighbor_mrgn/100 + 1;       % Scale values to (1 +- % margin)
				
				% Update SNR value - not here. not explicitly. This will be achieved implicitly through the
				% operation of adding more EMCs into the experimental vector. These additional EMCs will come
				% out of the use of more non-zero weights in the weights vector.
				% snr = snr * sqrt(N_neighbors+1);
				
			end
			

			% Store data in the current instance :
			%----------------------------------------------------------------------
			obj.experimental.nExprT2                   = length(loc);
			obj.experimental.neighbor_var			   = neighbor_var;
			%----------------------------------------------------------------------
			
		end
        function obj=simulate(obj,type)

            firstTE                         = obj.DBprop.firstTE;
            emcT2_full                      = obj.DBprop.emcT2_full;
			axis_T2_full                    = obj.DBprop.axis_T2_full;
			n_full_T2s                      = length(axis_T2_full);
            ETL                             = size(emcT2_full,ndims(emcT2_full));
						
			% IF TRUE -> Use the sparse EMC to build the fitting operator
			if strcmp(type,'forFit')
				
				d_T2                            = obj.operator.d_T2;
				maxT2                           = obj.operator.maxT2;
				minT2                           = obj.operator.minT2;
				sparsOrder                      = obj.operator.sparsOrder;
				min_T2_distance                 = obj.operator.min_T2_distance;

				n_sprs_T2s = abs(sparsOrder)*ETL;
				% add
				if (n_sprs_T2s) > n_full_T2s
					error('nT2s_to_fit_to is too large -- more than exists in EMC DB');
				end
				
				% Create a sparse T2 axis based on a given set of params
				T2axisSprs = calc_equispaced_T2(minT2,maxT2,firstTE,n_sprs_T2s,min_T2_distance,1,0);
				
				% Re-locate T2 values in the sparse axis to match the ones in the full axis
				T2axisSprs_matched = [];
				for idx = 1:length(T2axisSprs)
					[~,loc] = min(abs(round(T2axisSprs(idx)) - round(axis_T2_full)));
					if (loc > length(T2axisSprs_matched))
						T2axisSprs_matched(loc) = axis_T2_full(loc);
					end;
				end;
				T2axisSprs_matched = T2axisSprs_matched(T2axisSprs_matched~=0);
				T2axisSprs         = T2axisSprs_matched;
				axisT2             = T2axisSprs;
				n_T2s              = n_sprs_T2s;
				
				% Retrive sparsed EMC from the full EMC matrix
				extracted_locs = [];
				for idx = 1:length(axisT2)
					[min_val,lc] = min(abs(axisT2(idx)-axis_T2_full));
					if (min_val > min_T2_distance)
						error('Desired T2 weights vector contains a T2 location which does not exist in the simulated EMC DB T2 values');
					end
					extracted_locs = [extracted_locs lc];
				end
				emcT2_sprs = emcT2_full(:,:,:,extracted_locs,:);
				
				% [nbe] take following function code out to this function
				[E,~]=preprocess_emc(axisT2,emcT2_sprs,firstTE,[],[],ETL,type);
				
				%----------------------------------------------------------------------
				% Store operator data in the current instance :
				%----------------------------------------------------------------------
				obj.operator.axisT2            = axisT2;
				obj.operator.n_T2s             = n_T2s;
				obj.operator.emc               = emcT2_sprs;
				obj.operator.d_T2              = d_T2;
				obj.operator.minT2             = minT2;
				obj.operator.maxT2             = maxT2;
				obj.operator.sparsOrder        = sparsOrder;
				obj.operator.min_T2_distance   = min_T2_distance;
				obj.operator.E                 = E;
				%----------------------------------------------------------------------
				
				% IF FALSE -> Use all EMCs as an operator for fitting.
			elseif strcmp(type,'exprmt')
				snr                         = obj.experimental.snr;
				denoise_flag.denoise_op     = obj.experimental.denoise_flag.denoise_op;
				denoise_flag.r              = obj.experimental.denoise_flag.r;
				axisT2						= axis_T2_full;
				n_T2s					    = size(emcT2_full,4); % verify with Noam n_full_T2s ?
				emcT2_sprs					= emcT2_full;
				
				[~,E_noisy]=preprocess_emc(axisT2,emcT2_sprs,firstTE,snr,denoise_flag,ETL,type);
				
				%----------------------------------------------------------------------
				% Store experimental data in the current instance :
				%----------------------------------------------------------------------
				obj.experimental.E_noisy    = E_noisy;
				obj.experimental.n_T2s      = n_T2s;
				%----------------------------------------------------------------------
			end;	
			
		end
        function obj=genrtWeights(obj)
            
			axis_T2_full		= obj.DBprop.axis_T2_full;
			n_full_T2s			= length(axis_T2_full);
			relativeFract		= obj.experimental.relativeFract;
			loc					= obj.experimental.loc;
            winWidth			= obj.experimental.winWidth;
			mixMode				= obj.experimental.mixMode;
			N_neighbors         = obj.experimental.n_Neighbors;
			
            % ----------------------------------------------------------------------------------------
            % Prepare a multi T2 "experimental" vector from the noised EMC (E_noisy)
            % ----------------------------------------------------------------------------------------
            if (sum(round(relativeFract*1e5))~=1e5),  error('Sum fractions ~= 1');  end;

            weights = zeros(1,n_full_T2s);

			for neighbor_idx = 0:N_neighbors                 % Start at '0' because 0 pixel is the center one
				for loc_idx     = 1:length(loc)
					cur_T2_val  = loc(loc_idx);
					cur_fract   = relativeFract(loc_idx);
					
					if (neighbor_idx > 0)
						cur_T2_val = round(cur_T2_val * (neighbor_idx));
						cur_fract  =       cur_fract  * (neighbor_idx);
					end;
					% 		fprintf('cur_T2_val, fraction = %3.3f,  %3.3f\n',cur_T2_val, cur_fract);
					switch mixMode
						case 1           % Series of single T2 values
							[~, tmp_loc] = min(abs(cur_T2_val - round(axis_T2_full)));
							weights(tmp_loc) = weights(tmp_loc) + cur_fract;
							
						case 2           % Series of gaussians distribution around each T2
							T2low_tmp   = cur_T2_val - ceil(winWidth(loc_idx)/2);
							T2high_tmp  = cur_T2_val + ceil(winWidth(loc_idx)/2);
							[~, loc_low ] = min(abs(T2low_tmp -round(axis_T2_full)));
							[~, loc_high] = min(abs(T2high_tmp-round(axis_T2_full)));
							while ((loc_high-loc_low) < 3) && (loc_high<size(weights,2)) && (loc_low>=0)
								loc_low  = loc_low -1;
								loc_high = loc_high+1;
							end;
							T2_pk = transpose(gausswin((loc_high-loc_low+1),3));
							T2_pk = T2_pk / sum(T2_pk);
							weights(loc_low:loc_high) = weights(loc_low:loc_high) + T2_pk * cur_fract;
					end
				end
				
				% Save the "clean", center pixel, set of weights
				if (neighbor_idx == 0)
					weights_cntr_px = weights;
				end
				
				
			end
            
			%----------------------------------------------------------------------
			% Store data in the current instance : 
			%----------------------------------------------------------------------
			
            obj.experimental.weights         = transpose(weights         / sum(weights        ));
            obj.experimental.weights_cntr_px = transpose(weights_cntr_px / sum(weights_cntr_px));
            %----------------------------------------------------------------------
			
        end; 
        function obj=genrtExpEMC(obj)
			
			W			= obj.experimental.weights;
			W_cntrPix   = obj.experimental.weights_cntr_px;
			E_noisy     = obj.experimental.E_noisy;
			
            experimental_emc           = (squeeze(E_noisy))' * W;
            experimental_emc_center_px = (squeeze(E_noisy))' * W_cntrPix;
            
            % The experimental vector is generated from a noisy EMC DB
            obj.experimental.e         = experimental_emc;
            obj.experimental.e_cntrPix = experimental_emc_center_px;
			
		end

	end
end

