
function [Dictionary]=generateEMCdictionary(DBname)
% ____________________________________________________________________________
% DESCRIPTION : 
% This function creates an EMC dictionary according to a given set of
% parameters and a selected database file
% INPUT     - A path for the database '.mat' file (can be manually selected)
% OUTPUT    - The dictionary is saved in a '.mat' format to a selected folder
% ____________________________________________________________________________

%% Set parameters
p.minT2					= 1;	% [ms]
p.maxT2					= 5e2;	% [ms]
p.nT2					= 70;	% [#]
p.min_T2_dist			= 1;	% [ms]
p.T2base				= 20;	% [ms]
p.round_f				= 1;

p.minHeight				= 0;	% [#]
p.maxHeight				= 1;	% [#]
p.dH					= .1;	% [#]

p.nMaxComp				= 3;	% Number of T2 components in each dictionary item
p.plotFlag				= 0;

isGenerateWeights	    = 0;    % '1' for genertaing a dictionary / '0' for loading an existing one 


%% Load EMC database
orgCD       = cd;
compName    = getComputerName();
setCD(compName);
%cd(orgCD);

%Optional DB file 'C:\Users\me\Google Drive\Noam Omer\MultiCompT2\Data\2017_06_15_TAU_reference_32Echoes_400um_L_20mm_ABS.mat';
if (~exist('DBname','var'))
    cd ../
    cd ./DB
    [DBname,DBfolder]       = uigetfile(pwd,'Select a database file (.mat'')'); 
    disp('Database was loaded');
end

simData					= load(fullfile(DBfolder,DBname));             % Nz=4000, EMC_Abs_flag=1, TE=12, ETL=40, RefAng=120, BW=200
T2axisfull				= 1e3*simData.T2_tse_arr;	% converison from Sec to mSec
%emcT2B1full				= abs(simData.echo_train_modulation);
emcT2B1full				= abs(simData.echo_train_modulation_img);
firstTE					= 1e3*simData.TE_arr(1);	% converison from  Sec to mSec
TE_arr					= simData.TE_arr;
B1_1					= find(simData.B1_scaling_arr==1); % curently implemented to constant B1
emcT2full               = emcT2B1full(1,B1_1,:,:,:);
T2axisLog               = calc_equispaced_T2(p.minT2,p.maxT2,firstTE,p.nT2,p.min_T2_dist,p.round_f,p.plotFlag);
p.T2axis                = T2axisLog;
p.B1					= B1_1;


%% Weights permutations matrix

% Generate a new weights dictionary
if isGenerateWeights
    disp('Generating weights permutations matrix ...')
    tic
    [weights] = generateWeights(p);
    time2createDictionary=toc;
    disp(sprintf('Time to generate the weights permutations matrix : %.2f minutes',time2createDictionary./60))

% Load an existing weights permutation matrix from directory
else
    orgCD=cd;
    cd ../
    cd ./Results/WeightsPermutations
	file2load=uigetfile();
    load(file2load);
    cd(orgCD);
end

% Scale Dictionary to 0..1
 weights=weights./p.maxHeight;

%% Generate the dictionary

% Retrive sparsed EMC from the full EMC matrix
extracted_locs = [];
for idx = 1:length(T2axisLog)
    [min_val,lc] = min(abs(T2axisLog(idx)-T2axisfull));
    if (min_val > T2axisLog)
        error('Desired T2 weights vector contains a T2 location which does not exist in the simulated EMC DB T2 values');
    end
    extracted_locs = [extracted_locs lc];
end


emcT2sparse = emcT2full(:,:,:,extracted_locs,:);
if length(size(emcT2full))==4 % A single T1 database
	emcT2sparse = normalize_4D_EMC(emcT2sparse);
else 
	emcT2sparse = normalize_5D_EMC(emcT2sparse);
end

[n_B0,n_B1,n_T1,n_T2,ETL]	= size(emcT2sparse);

% Calculate a decay matrix for each simulated EMC vector
% accorting to the first TE (instead of the EMC starting at '1')
decayVec         = transpose(exp(-firstTE./(1e3*T2axisLog)));
decayMat         = repmat(decayVec,1,ETL);
emcT2sparseDecay = squeeze(emcT2sparse) .* decayMat;

if size(weights,1)~=size(emcT2sparseDecay,1)
	error('The number of T2 values within the selected weights matrix do not match to DB')
end
emc_curves = transpose(emcT2sparseDecay)*weights; % size: [nT2 X nWeights]

%% Store data 

Dictionary.T2axis		= T2axisLog;
Dictionary.TE_arr		= TE_arr;
Dictionary.weights		= weights;
Dictionary.curves		= emc_curves;
Dictionary.params		= p;
Dictionary.DPpath       = DBname; 
Dictionary.emcNormDec   = emcT2sparseDecay;

disp('Dictionary was created succssesfuly !');

% Save dictionary EMC curves
orgCD=cd;
cd ../
cd ./Results/Dictionary
strr=datestr(datetime);
strr([strfind(strr,' '),strfind(strr,':')])='_';
save(['Dictionary',strr],'Dictionary');
cd(orgCD);
disp(['Dictionary saved under: ','Dictionary',strr]);
return;




 