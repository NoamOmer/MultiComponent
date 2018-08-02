
clcl;clear all;

% Load Dictionary

%% Main Simulator Script
clc;clcl;

% Initalizaiton
%----------------------------------------------------------------------
% Set current directory
compName=getComputerName;
setCD(compName);
cd ..
cd Results\Dictionary	

params = setParams();		

orgCD=cd;
cd ../
cd ./Results/Dictionary
file2load=uigetfile();
load(file2load);
cd(orgCD);
curves=transpose(Dictionary.curves);
nCurves=size(curves,1);


% Load MSE data
%----------------------------------------------------------------------
cd ../../../
%cd Scans\SEMC149
cd Scans\27_12_17-Dvir\MSE
SEMCDataDir=cd;
vendor='Siemens';
opts.PD          = 1;
opts.T2          = 1;
opts.dcm         = 1;
opts.nii         = 0;
opts.Siemens_dir = [];

switch vendor
	case 'Siemens'
		[dcm_info,dcm_info_parsed,im_SEMC_mSl,slices_info] = Read_mSl_Data(SEMCDataDir);
	case 'GE'
		[dcm_info,dcm_info_parsed,im_SEMC_mSl,slices_info] = Read_mSl_Data_GE(SEMCDataDir);
end

if  0%(Dictionary.TE_arr(1) ~= dcm_info_parsed.alTE(2)*1e-6) [NO]
	error('dictionary do not compatible to the data -  TE mismatch');
else
	firstTE  = dcm_info_parsed.alTE(2)*1e-6; % [sec]
end


% Normalize the measured signals
sigMat=squeeze(im_SEMC_mSl);
normSigs = normalize_SEMC(sigMat);

% From this stage the process is done for the normalize signals 
sigs=normSigs; 

% Get the desired ROI - currently manually 
%(in the future this stage will be performed based on segmentation)
numSlice2plot=floor(size(sigs,3)/2);
[~,limits] = imcrop(sigs(:,:,numSlice2plot));
limits=round(limits);
close;

% Extract only the signals from a given segmet segment 
indExtactedRows=limits(1)+(1:limits(3));
indExtactedCols=limits(2)+(1:limits(4));

sigPerSegment=zeros(size(indExtactedRows,2),size(indExtactedCols,2),size(sigs,3));
for i_slice=1:size(sigs,3)
	sigPerSegment(:,:,i_slice) = sigs(indExtactedRows,indExtactedCols,i_slice);
end
sigPerSegment=reshape(sigPerSegment,size(sigPerSegment,1)*size(sigPerSegment,2),[]);
nSigs=size(sigPerSegment,1);
% mask=zeros(size(sigs(:,:,1)));
% mask(intMat)=1;
%mask_SEMC = CreateSegmentMask(sigs(:,:,1),intMat);



% Calculate scores
%----------------------------------------------------------------------

nEchoes4Score=10;
echosVec=1:nEchoes4Score;


scores=zeros(nSigs,size(curves,1));
h_wb = waitbar(0,'Calculating scores...');

for iSig = 1:size(sigPerSegment,1)
	tmpSig=sigPerSegment(iSig,:);
	tmpSigMap=repmat(tmpSig,size(curves,1),1);
	scores(iSig,:)=sqrt(sum(power(tmpSigMap(:,echosVec)-curves(:,echosVec),2),2)); 
	waitbar(iSig / nSigs,h_wb,sprintf('processing signals ... %d / %d',iSig,nSigs));
end
close(h_wb);

sortedScores=sort(scores,'descend');
scoreThr=0.01;
candidates=nan(nCurves,2);

% OPTION 2
vxcvcx+1;
nBasisElements=3;
basisElementsInd=nSigs.*ones(nBasisElements,3);
counter=0;
for iCurve=1:nCurves
    [loc]=find(sortedScores(:,iCurve)<=scoreThr,1,'first'); % Change this so the maximal value will be replaced each time
    if ~isempty(loc) && sum(basisElementsInd(:,3)>loc)
        if counter==nBasisElements
            replacementInd=find(basisElementsInd(:,3)>sortedScores(loc,iCurve),1,'first');
            
            if ~isempty(replacementInd)
                basisElementsInd(replacementInd,1:3)=[iCurve,sortedScores(loc,iCurve),loc];
            end
            
        else
            counter=counter+1;
            basisElementsInd(counter,1:3)=[iCurve,sortedScores(loc,iCurve),loc];
        end

    end
end




% Plot some sorted scores for example 
figure();semilogy(sortedScores(:,[1:5e3:end]),'.-');


% OPTION 1 
% for iCurve=1:nCurves
%     [loc]=find(sortedScores(:,iCurve)<=scoreThr,1,'first');
%     if loc
%         candidates(iCurve,1)=loc;
%         candidates(iCurve,2)=sortedScores(loc,iCurve);
%     end
% end



% figure();
% for ii=(1e5+(1:1000))
% 	plot(1:length(tmpSig),tmpSig,1:20,curves(ii,1:20),'.');
% 	hold on;
% end


% [~ ,curveInd]=nanmin(candidates);
% tmpCandidates=candidates;tmpCandidates(curveInd)=NaN;
% [~ ,curveInd]=nanmin(tmpCandidates);



% Normalize Curve!!!!! Bug !
% 
% xmin = 0; % lower bound
% xmax = nSigs; % upper bound
% scoresResolution=10;
% h = nSigs/scoresResolution; % bin width
% edges = xmin:h:xmax;
% 
% % Allocate scores in fixed width bins 
% [counts,centers] = histc(scores(:,i_curve),edges);
% counts=counts./nSigs;
% 
% [maxVals,maxInds]=max(counts);
% [maxVals,maxInds]=max(counts);
% 
% figure();
% for i=500:550
%     subplot(1,2,1);semilogy(counts(:,i)','.-');hold on;
%     subplot(1,2,2);plot(counts(:,i)','.-');hold on;
% end
% 
% 
% 
% basisElements=find(sum(counts(1:2,:))>.984848478);
% figure();plot(basisElements);
% 
% 
% 
% 
% 
