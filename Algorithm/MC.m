function main()

clcl;
 path='C:\Users\me\Google Drive\Noam Omer\MultiCompT2\2016_09_22_T2_1_1000_B1_100_for_MultiT2Fit_400um_L_20mm_ABS.mat';
% simData						= load(path);  
% axisT2						= 1e3*simData.T2_tse_arr; % converison from mSec to Sec
% emcT2						= abs(simData.echo_train_modulation);
% [nB0,nB1,nT1,nT2,ETL]		= size(emcT2);
% firstTE						= 1e3*simData.TE_arr(1);	% converison from mSec to Sec1
% TE_arr						= simData.TE_arr;
% E = squeeze(emcT2(1,1,:,:,:));

load('G:\My Drive\Noam Omer\MultiCompT2\mat\e.mat');
load('G:\My Drive\Noam Omer\MultiCompT2\mat\E.mat');

% nETL	     = min(size(E));
% firstTE    = 10;		%	  [ms]
% axisT2     = 1:nT2;		% T2s [ms]
% decayVec   = transpose(exp(-firstTE./(axisT2)));
% decayMat   = repmat(decayVec,1,nETL);

nT2					= max(size(E));
nETL				= min(size(E));
dw				    = 0.01;
nExprmnt		    = 1e3;
nVox			    = 1; % currently one signal is tested
plotFlag		    = 1;
similarityTestType  = 2; 
randomnessType		= 1;

% Allocate 
scores	= zeros(nExprmnt,nVox);
w		= zeros(nT2,nExprmnt);
w_dig	= zeros(size(w));
er		= zeros(nETL,nExprmnt);




%% Randomness
rng('shuffle')
switch randomnessType
	case 1 % option 1 :
		tmp = rand(nT2,nExprmnt);		% Uniformly distributed random numbers
		w = bsxfun(@rdivide,tmp,sum(tmp,2));
	case 2 % option 2 :
		w = randfixedsum(nT2,nExprmnt,1,0,1);
end


%%

%C = bsxfun(@digitizeWeights,w,dw);
for i_exprmnt=1:nExprmnt
	w_dig(:,i_exprmnt)=digitizeWeights(w(:,i_exprmnt),dw);
	er(:,i_exprmnt)=E'*w_dig(:,i_exprmnt);
	scores(i_exprmnt,1) = calcScore(er(:,i_exprmnt),e,similarityTestType);
end

[~,ind]=max(scores(:,1));
figure(); plot(1:length(er(:,ind)),er(:,ind),1:length(e),e);grid on;


end



function score = calcScore(exprmtEMC,EMC,type)
	switch type 
		case 1  % Correlation
			score = corrcoef(exprmtEMC,EMC);
			score = score(1,2);
		case 2	% Dynamic Time Warping
			score=dtw(exprmtEMC,EMC);
		case 3
	end
end



function w = digitizeWeights(w,dw)

nBins=(1/dw);						% Number of weights bins
wBins=linspace(dw/2,1-dw/2,nBins);  % Get number of bins

for i=1:length(w)
	
	if w(i)< wBins(1)
		w(i) =0;
		continue;
	end
	
	if w(i)>= wBins(end)
		w(i) =1;
		continue;
	end
	
	for j=1:length(wBins)
		if (w(i)>=wBins(j) && w(i)<wBins(j+1))
			w(i)=mean([wBins(j),wBins(j+1)]);
		end
	end
	
end

end