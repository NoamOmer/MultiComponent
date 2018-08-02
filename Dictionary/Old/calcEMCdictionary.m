
function [Dictionary]=generateEMCdictionary()


% Set parameters

p.minT2        = 1;	% [ms]
p.maxT2        = 1e3;	% [ms]
p.dT2          = 1;	% [ms]
p.minHeight    = 0;	% [#]
p.maxHeight    = 10;	% [#]
p.dH           = 1;	% [#]
p.nT2          = 50;	% [#]
p.min_T2_dist  = 2;	% [ms]
p.T2base       = 40;	% [ms]
p.round_f      = 1;
p.nMaxComp     = 3;	% Number of components within dictionary
p.plotFlag     = 0;
isGenerateDictionary= 0;     % '1' for genertaing a dictionary / '0' for loading an existing one 


% Load EMC database
simData					= load('2016_09_22_T2_1_1000_B1_100_for_MultiT2Fit_400um_L_20mm_ABS.mat');      % Nz=4000, EMC_Abs_flag=1, TE=12, ETL=40, RefAng=120, BW=200
axisT2full				= 1e3*simData.T2_tse_arr;	% converison from Sec to mSec
emcT2full				= abs(simData.echo_train_modulation);
firstTE					= 1e3*simData.TE_arr(1);	% converison from mSec to Sec
TE_arr					= simData.TE_arr;
B1_1					= find(simData.B1_scaling_arr==1);
emcT2full               = emcT2full(1,B1_1,:,:,:);
T2axis = calc_equispaced_T2(p.minT2,p.maxT2,firstTE,p.nT2,p.min_T2_dist,p.round_f,p.plotFlag);


% Generate a new dictionary
if isGenerateDictionary
    tic
    [weights,T2axis] = generateWeights(params);
    time2createDictionary=toc;

% Load Dictionay
else
    orgCD=cd;
    cd ../
    cd ./results/Dictionary
    load('weights.mat');
    cd(orgCD);
end

% Normelaize Dictionary
 weights=weights./p.maxHeight;

	
% Re-locate T2 values in the sparse axis to match the ones in the full axis
T2axisSprs_matched = [];
for idx = 1:length(T2axis)
    [~,loc] = min(abs(round(T2axis(idx)) - round(axisT2full)));
    if (loc > length(T2axisSprs_matched))
        T2axisSprs_matched(loc) = axisT2full(loc);
    end
end
T2axisSprs_matched = T2axisSprs_matched(T2axisSprs_matched~=0);
T2axis         = T2axisSprs_matched;


% Retrive sparsed EMC from the full EMC matrix
extracted_locs = [];
for idx = 1:length(T2axis)
    [min_val,lc] = min(abs(T2axis(idx)-axisT2full));
    if (min_val > T2axis)
        error('Desired T2 weights vector contains a T2 location which does not exist in the simulated EMC DB T2 values');
    end
    extracted_locs = [extracted_locs lc];
end

emc = emcT2full(:,:,:,extracted_locs,:);
emcNorm = normalize_5D_EMC(emc);

[n_B0,n_B1,n_T1,n_T2,ETL]	= size(emc);
% Calculate a decay matrix for each simulated EMC vector
% accorting to the first TE (instead of the EMC starting at '1')
decayVec        = transpose(exp(-firstTE./(1e3*T2axis)));
decayMat        = repmat(decayVec,1,ETL);
emcNormDecay = squeeze(emcNorm).*decayMat;

% Reshape to concatenate all EMC cerves into a 2D matrix [ETL x nT2]
E = permute(emcNormDecay,[5 4 2 3 1]);
E = reshape(E,40,[],1);

curves = E*weights;

% Store data in the output struct 
Dictionary.T2axis=T2axis;
Dictionary.TE_arr=TE_arr;
Dictionary.weights=weights;
Dictionary.curves=curves;
Dictionary.params=p;

% Save dictionary EMC curves
orgCD=cd;
cd ../
cd ./results/Dictionary
strr=datestr(datetime);
strr(strfind(strr,' '))='_';
strr(strfind(strr,':'))='_';
save(['Dictionary',strr],'Dictionary');
cd(orgCD);
end



