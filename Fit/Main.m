%% Main Simulator Script
clc;clcl;

% Initalizaiton
%----------------------------------------------------------------------
% Set current directory
compName=getComputerName;
setCD(compName);
params = setParams();			

% Create fitting operator -> E
%----------------------------------------------------------------------
type= 'forFit';
fitOpr = Simulator;                     % create a new instance of type Simulator
fitOpr = fitOpr.getParams(type);		% get parameters for synthesizing a fitting operator
fitOpr = fitOpr.loadSimData();          % load acquisition operator (E_full)
fitOpr = fitOpr.simulate(type);			% generate fitting operator (E)
	

% Create an experimental EMC -> e
%----------------------------------------------------------------------
type   = 'exprmt';

%for ii=1:6
	
expmnt = Simulator;                   % create a new instance of type Simulator
expmnt = expmnt.getParams(type);	  % get parameters for synthesizing experimental EMC
% [nbe] don't re-load sim data - no need
expmnt = expmnt.loadSimData();        % load acquisition operator (E_full)
expmnt = expmnt.simulate(type);       % generate a noisy experimental operator (E_noisy)
expmnt = expmnt.genrtWeights();		  % generate vector of weights for the experimental EMC  
expmnt = expmnt.genrtExpEMC();		  % generate experimental emc 	
% [NBE] the last 3 functions should be combined into one single funciton

%eval(['expmnt',num2str(ii),'=expmnt;'])

%end
% Solver - Reconstructing the unknown weights vector (w) 
%----------------------------------------------------------------------
E      = fitOpr.operator.E';		% Fitting operator (E)
e      = expmnt.experimental.e;		% Experimental vector (e)


% Solve !
%----------------------------------------------------------------------
sol    = solver(E,e,params);		 


% Plot Reconstruction Results
% ----------------------------------------------------
emcFittedNorm	= normalizeEMC(squeeze(fitOpr.operator.emc));
eFitted			= transpose(emcFittedNorm)*sol;

Displayer.plotEMCcomparison(expmnt,eFitted)
Displayer.plotWeightsComparison(fitOpr,expmnt,sol);

%%





