%% Main Algorithm Function

%clc;
clcl;
go2myRootDir

% Set current directory
%----------------------------------------------------------------------
compName=getComputerName;
setCD(compName);
rootDir=cd;


% Configuration
%----------------------------------------------------------------------
nEchoes4Cost	= 15;               % Number of echoes used for cost's calculation ([1:nETL])
precentThr		= .05;              % A maximal relative cost cutoff for a T2-weights-vector to be a basis element candidate
nBasisElements	= 15;				% Number of basis elements to be included in the analysis
findCostThr     = 2;                % Type of method to find the cost threshold: '1' Sort all weights of all signals \ '2' Histogram based detection
echosVec		= 1:nEchoes4Cost;
sigLevel4cost   = .1;               % A threshold for the minimal signal level to be included in the cost calculation


% Solver Parameters
%----------------------------------------------------------------------
solType                                  = 'singelVoxel';%'singelVoxel'; %'singelVoxel';   \ 'E'
fitParams.fit_TYPE                       = 2;     % [1: Linear  ||  2: Quadratic  || 3: OMP]
fitParams.fit_nRep                       = 1;     % Repeat fit in case it can use prior info from previous runs
fitParams.Sum_fitted_w_Eql_to_one_flag   = 1;

fitParams.Tikh_reg_flag                  = 0;     % [0 | 1] Enable / disable the use of Tikhonov regularization in fitting
fitParams.lambda_Tikh                    = 1e-6;  % Tikhonov regularization weight
fitParams.w0_Tikh                        = 0;     % Apply Tikhonov around non zero baseline value

fitParams.lambda_L1                      = 1e+7;  % L1 regularization weight

% NOTES:
% quadprog with nRep=2 seems to de-noise the vector of weights
%fitParams.neighbor_mrgn                  = 15;    % [%] variability of T2 value and relative fraction between neighboring pixels ~5%
%fitParams.N_neighbors                    = 0;     % Number of neighboring pixels (0: no neighbors | 8: 3x3 kernel | 15: 4x4 kernel)


% Load EMC dictionary
%----------------------------------------------------------------------
cd ..
cd Results\Dictionary

scanner = menu('Data from which sccanner should be analyzed ?','Bruker','Siemens');

switch scanner
	case  1
		cd Bruker
		file2load	= uigetfile('','Load a dictionary file...'); 
		load(file2load); cd(rootDir);
	case 2
		cd Siemens
		file2load	= uigetfile('','Load a dictionary file...'); 
		load(file2load); cd(rootDir);
end



curves		= transpose(Dictionary.curves);     % [nPermutations X nETL]
weights		= transpose(Dictionary.weights);    % [nPermutations X nT2]
nCurves		= size(curves,1);                   % [nPermutations]
T2axis		= Dictionary.T2axis;                % [nT2]
TEarr		= Dictionary.TE_arr;				% [nETL]
nT2			= length(T2axis);


%% Analysis Modes
dataType = 'realData';%'simulatedData';%'simulatedData';%;%; %

d = dataParser(dataType);

nSegments	= d.nSegments;
nSigsPerSeg = d.nSigsPerSeg;
indVoxPerSeg= d.indVoxPerSeg;
normSigs	= d.normSigs;
sigs		= d.sigs;
nRows		= d.nRows;
nCols		= d.nCols;
nSlice		= d.nSlice;
nETL		= d.nETL;
firstTE		= d.firstTE;

% Verify that a correct analysis is performed
if nETL~=length(TEarr)
	error('The number of echoes in dictionary (TEarr) and the data (nETL) missmatch !');
end

if firstTE~=TEarr(1) 
	error('The TE in dictionary and TE of the simulated data mismatch !');
end


% 4) Analyze each segment
for iSeg = 1:nSegments
	
	% Initiate for the current segment
	nSigsPerSeg			= size(indVoxPerSeg{iSeg},1);			% Number of voxels in a segment
	locInd				= cell(nSigsPerSeg,2);					% Allocate a cell for the segment voxels spatial indexes (1st column) and its T2 weights vector (2nd Column)
	for i_Vox=1:nSigsPerSeg
		locInd{i_Vox,1}		= indVoxPerSeg{1,iSeg}(i_Vox);		% Store the spatial indexes of all the segment voxels
	end
	sigPerSegment		= zeros(nSigsPerSeg,nETL);				% Allocate the signals of the current segment
	basisElementsInd	= nSigsPerSeg.*ones(nBasisElements,3);	% Initiate the matrix indicating the basis elements indexes with the highest index
	
	% Extract only the signals of a given segment
	sigPerSegment		= normSigs(indVoxPerSeg{iSeg});
	
% 	for iSlice=1:nSlice
% 		sigPerSlice = sigs(:,:,iSlice);
% 		sigPerSegment(:,iSlice) = sigPerSlice(indVoxPerSeg{iSeg});
% 	end
	
	% Calculate costs
	%----------------------------------------------------------------------
	try
		costs	= zeros(nSigsPerSeg,nCurves);			% Allocate a costs matrix (nVoxels X nWeights)
	catch
		costs	= zeros(nSigsPerSeg/2,nCurves);
	end
	
	h_wb	= waitbar(0,'Calculating costs...');	% Open a waitbar display
	
	% Calculte the cost based on L2 norm of differences
	
	
	
	for iSig = 1:nSigsPerSeg
		tmpSig			= sigPerSegment(iSig,:);
		tmpSigMat		= repmat(tmpSig,size(curves,1),1);
		% Option 1 - direct calculation of the root of the sum of square differences
		% Note : Do not include signals under a threshold 'sigLevel4cost' into the cost calculation
		% 				tmpSigMap(tmpSigMap<=sigLevel4cost)= NaN;
		
		tmpSigMat(tmpSigMat<=sigLevel4cost)= 0;
		curves   (curves   <=sigLevel4cost)= 0;
		curves(tmpSigMat==0)=0;
		tmpSigMat(curves==0)=0;
		costs(iSig,:)	= sqrt(sum(power(tmpSigMat(:,echosVec)-curves(:,echosVec),2),2));
		
		
		cw                  = costs(iSig,:);
		
		
		% Find the weight with the minimal score
		[can,indCan]	= min(cw);
		can				= weights(indCan,:);
		disp('Minimal cost is not 0');
		
		if strfind(dataType,'simulatedData')
			
			% Find the weight with score '0'
			desiredWeightInd    = getDesiredWeight(weights,sim.experimental.loc,sim.experimental.relativeFract,T2axis);
			can                 = weights(find(cw<2*eps),:);

			% Find the weight with the minimal score
			if ~isempty(can)
				indCan=desiredWeightInd;
			end
			% In a single simulated data case  there is no need to calculate the costs
			costs =repmat(costs(iSig,:)  ,nSigsPerSeg,1)  ;
			break;
		end
		
		% Option 2 - using 'norm' funciton
		%norms          = cellfun(@norm,num2cell((tmpSigMap(:,echosVec)-curves(:,echosVec)),2));
		% Option 3 -Use vecnorm starting from ver17b
		%costs(iSig,:)  = vecnorm((tmpSigMap(:,echosVec)-curves(:,echosVec))',2);
		
		waitbar(iSig / nSigsPerSeg,h_wb,sprintf('processing signals ... %d / %d',iSig,nSigsPerSeg));
	end
	close(h_wb);
	
	% Figure 1 : plot the weight with the lowest cost
% 	figure(1);
% 	subplot(2,1,1);semilogx(repmat(T2axis,size(can,1),1),can,'--');grid;
% 	title('Weight with the lowest cost');
% 	hold on;subplot(2,1,1);semilogx(T2axis,weights(desiredWeightInd,:));
% 	legend('Candidate','Ground Truth');
% 	tmpSig(tmpSig   <=sigLevel4cost)= 0;
% 	
% 	subplot(2,1,2);plot(TEarr,tmpSig,'--');grid;xlabel('ms');
% 	hold on; subplot(2,1,2);
% 	plot(TEarr,sim.experimental.experimental_denoised_emc,'.-');grid on;xlabel('ms');
% 	legend('Candidate','Ground Truth');
% 	disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~');
% 	disp('Last trail results:');
% 	disp(['SNR = ',num2str(sim.experimental.snr)]);
% 	disp('Experimental : ');
% 	disp([sim.experimental.loc,sim.experimental.relativeFract]);
% 	disp('Minimal Score Weight : ');
% 	disp([T2axis(find(can)),weights(indCan,find(weights(indCan,:)))])
% 	disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~');
	
	% 			normCosts= (max(costs(:))-costs)./(max(costs(:))-min(costs(:)));
	% 			sum(normCosts.^2);
	% 			indSpatialProximity=[1,2,3,4,5,6];
	% 			normCostsProximity=sum(normCosts(indSpatialProximity,:));
	% Find the cost threshold based on all the costs within the the segment
	switch findCostThr
		
		case 1
			% Option 1 : Sort all costs of all signals - accurate but might take a lot of time !
			precentAxes				= (1:nCurves*nSigsPerSeg)/(nCurves*nSigsPerSeg);
			allCostsSortedVec		= sort(costs(:),'ascend');
			accumulatedCostsVec	= cumsum(allCostsSortedVec);
			[minVal,indMinVal]		= min(abs(repmat(precentThr,1,(nCurves*nSigsPerSeg))-(precentAxes)));
			costThr				= allCostsSortedVec(indMinVal);
			
		case 2
			% Option 2 : Sort the histogram of the costs of all signals - less accurate but much faster !
			[N,edges]				= histcounts(costs(:));
			N_norm					= N./sum(N);
			accumulatedCostsVec		= cumsum(N_norm);
			[minVal,indMinVal]		= min(abs(repmat(precentThr,1,(length(N)))-(accumulatedCostsVec)));
			costThr					= edges(indMinVal); %(edges(indMinVal-1)+edges(indMinVal))/2;
			
	end
	
	% Search for basis elements
	sortedCosts = sort(costs,'descend'); % Sort the costs of each signal
	%figure;plot(sortedCosts(:,[1,1e3,2e3]));grid on;xlabel('Voxel');ylabel('Cost');legend('W1','W2','W3');
	modScore=zeros(1,size(sortedCosts,2));
	for iCurve=1:nCurves
		
		% Option 1: For each curve find the first index with a cost lower\equal to the selected threshold
		[loc]=find(sortedCosts(:,iCurve)<=costThr,1,'first');
		if ~isempty(loc)
			val=sortedCosts(loc,iCurve);
			modScore(iCurve) = val./((loc).^2);
		else
			modScore(iCurve)=NaN;
		end
		
		% Option 2: Keep the current weight if it has reached s cost which is less then the 'costThr'
		%			and in case it reaches this cost before the current basis elements candidates
		% if ~isempty(loc) &&  sum(basisElementsInd(:,2)>val)
		% if sum(basisElementsInd(:,3)-repmat(loc,nBasisElements,1))
		%[~,replacementInd]=max(basisElementsInd(:,3)-repmat(loc,nBasisElements,1));
		%basisElementsInd(replacementInd,1:3)=[iCurve,sortedCosts(loc,iCurve),loc];
		% else
		% [~,replacementInd]=max(basisElementsInd(:,2));
		% basisElementsInd(replacementInd,1:3)=[iCurve,sortedCosts(loc,iCurve),loc];
		% end
		% end
		
	end
	
	[modScoreSorted,modScoreInd]=sort(modScore);
	basisElementsInd=[modScoreInd(1:nBasisElements)',modScore(modScoreInd(1:nBasisElements))'];
	basisElements	= transpose(weights(basisElementsInd(:,1),:));
	E				= transpose(curves(basisElementsInd(:,1),echosVec));	% EMC Dictionary
	E_concate		= repmat(E,nSigsPerSeg,1);
	WV				= zeros(nT2,nSigsPerSeg);								% Unknown T2-weights matrix
	
	% Related to option 2:
	% 			basisElements	= transpose(weights(basisElementsInd(:,1),:));			% T2 basis elemnts weights
	% 			E				= transpose(curves(basisElementsInd(:,1),echosVec));	% EMC Dictionary
	% 			E_concate		= repmat(E,nSigsPerSeg,1);
	% 			WV				= zeros(nT2,nSigsPerSeg);								% Unknown T2-weights matrix
	
	% Solve the optimization problem for each voxel within the current segment
	for i_sig=1:nSigsPerSeg
		
		
		switch solType
			
			% A single voxel solver
			%----------------------------------------------------------------------
			case 'singelVoxel'
				e       = sigPerSegment(i_sig,echosVec);
				sol     = solver(E,e,fitParams); % The minimizer for the problem : ||E^W^-e||^2
				
				% Contrast2Noise Neighbors(4) solver
				%----------------------------------------------------------------------
			case 'contrast2noise'
				[indRow,indCol]=ind2sub(size(im.org),indVoxPerSeg{1,iSeg}(i_sig));
				upVox       = find(indVoxPerSeg{1,iSeg}==sub2ind(size(im.org),indRow+1,indCol));
				downVox     = find(indVoxPerSeg{1,iSeg}==sub2ind(size(im.org),indRow-1,indCol));
				rightVox    = find(indVoxPerSeg{1,iSeg}==sub2ind(size(im.org),indRow,indCol+1));
				leftVox     = find(indVoxPerSeg{1,iSeg}==sub2ind(size(im.org),indRow,indCol-1));
				eTmp        = sigPerSegment(i_sig,echosVec);
				
				neibrCount=1;
				if ~isempty(upVox);     eTmp=eTmp+sigPerSegment(upVox,echosVec);    neibrCount=neibrCount+1;    end;
				if ~isempty(downVox);   eTmp=eTmp+sigPerSegment(downVox,echosVec);  neibrCount=neibrCount+1;    end;
				if ~isempty(rightVox);  eTmp=eTmp+sigPerSegment(rightVox,echosVec); neibrCount=neibrCount+1;    end;
				if ~isempty(leftVox);   eTmp=eTmp+sigPerSegment(leftVox,echosVec);  neibrCount=neibrCount+1;    end;
				
				eAvg    = eTmp./neibrCount;
				sol     = solver(E,eAvg,fitParams); % The minimizer for the problem : ||E^W^-e||^2
				
				
				% Neighbors(4) Signals solver
				%----------------------------------------------------------------------
			case '4Neighbors'
				[indRow,indCol]=ind2sub(size(im.org),indVoxPerSeg{1,iSeg}(i_sig));
				upVox       = find(indVoxPerSeg{1,iSeg}==sub2ind(size(im.org),indRow+1,indCol));
				downVox     = find(indVoxPerSeg{1,iSeg}==sub2ind(size(im.org),indRow-1,indCol));
				rightVox    = find(indVoxPerSeg{1,iSeg}==sub2ind(size(im.org),indRow,indCol+1));
				leftVox     = find(indVoxPerSeg{1,iSeg}==sub2ind(size(im.org),indRow,indCol-1));
				eTmp     = sigPerSegment(i_sig,echosVec);
				
				if ~isempty(upVox);     eTmp=[eTmp;sigPerSegment(upVox,echosVec)];        end;
				if ~isempty(downVox);   eTmp=[eTmp;sigPerSegment(downVox,echosVec)];      end;
				if ~isempty(rightVox);  eTmp=[eTmp;sigPerSegment(rightVox,echosVec)];     end;
				if ~isempty(leftVox);   eTmp=[eTmp;sigPerSegment(leftVox,echosVec)];      end;
				
				eConcat = reshape(transpose(eTmp),1,[]);
				sol     = solver(repmat(E,size(eTmp,1),1),eConcat,fitParams); % The minimizer for the problem : ||E^W^-e||^2
				
		end
		
		WV(:,i_sig)		= (basisElements*sol);  % Final weights vector for the voxel
		locInd(i_sig,2)	= {WV(:,i_sig)};        % Store each T2 weights vector next to its orginal spatial location
		
	end
	
	WMat{iSeg,1}= {WV};
	
	% Plot Analysis Resutls
	%----------------------------------------------------------------------
	
	% Figure 2 : Basis elements
	figure(2);
	subplot(1,2,1);
	weights2plot	= weights(basisElementsInd(:,1),:);
	h_BasisWeights  = semilogx(repmat(T2axis',1,nBasisElements),weights2plot','.-');
	
	title(sprintf('T_2 Distributions of each basis element (%d)',nBasisElements));
	xlabel('T_2 Relaxation Time [mSec]');
	ylabel('Normalized Component Amplitude');
	axis([T2axis(1),T2axis(end),0,1.1]);
	legendStr={};
	for i_element=1:nBasisElements
		legendStr= [legendStr {['Basis Element ', num2str(i_element)]}];
	end
	legend(legendStr);
	grid on;
	
	subplot(1,2,2);
	h_BasisElements=semilogy(repmat(TEarr',1,nBasisElements),curves(basisElementsInd(:,1),:)','.-');
	legend(legendStr);
	grid on;
	for i_plot=1:nBasisElements
		h_BasisElements(i_plot).Marker='o';
		h_BasisElements(i_plot).MarkerSize=6;
	end
	xlabel('Time [mSec]');
	ylabel('Normalized Amplitude');
	title('EMC');
	
	
	% Figure 2 : Sorted costs for examples
	%figure(2);
	%subplot(2,1,1);semilogy(sortedCosts(:,[1:5e3:end]),'.-');grid;title('Sorted costs of random weights');
	%subplot(2,1,2);semilogy(sortedCosts(:,basisElementsInd(:,1)),'.-');grid;title('Sorted costs of selected basis elements');
	
	
	% Figure 3 : T2 Distribution for all voxels within the segment
	figure(3);
	[X,Y] = meshgrid(T2axis,1:nSigsPerSeg); Z = WV';
	Z(Z==0)=NaN;
	surf(X,Y,Z); %contourf(X,Y,Z)
	%axis([T2axis(1) 200 T2axis(1) nSigsPerSeg])
	view(2);
	xlabel('T_2 Relaxation Time [mSec]');
	ylabel('Voxel [#]');
	title(sprintf('T_2  Distribution - within a Segment (%d Voxels)',nSigsPerSeg));
	colorbar;
	
	%k = waitforbuttonpress; % Wait for the user click to move forward to the next segment
	
end


figure(4);
for iSeg=1:nSegments
	segMat=cell2mat(WMat{iSeg,1});
	subplot(1,1,iSeg);semilogx(repmat(T2axis',1,size(segMat,2)),segMat,'-o');grid on;
	hold on;
	semilogx(T2axis,weights(desiredWeightInd,:));grid on;
	axis([T2axis(1) 5e3 0 1.1]);
	legend('Candidate','Ground Truth');
end

title(sprintf('T_2  Distribution - within a Segment (%d Voxels)',nSigsPerSeg));


% Figure 5 : plot the final T2 distribution for each voxel
figure(5);
for i_plot=1:nSigsPerSeg
	subplot(4,ceil(nSigsPerSeg/4),i_plot);
	stem(T2axis,segMat(:,i_plot)','linestyle','-');
	grid on;
	set(gca,'xscal','log');
	set(get(gca,'XAxis'),'Limits',[1 T2axis(end)]);
end



% Figure 6 : plot a T2 map for several T2 values
% Option 1 : For brain data
% T2range2plot=[20 50 80 110 140 170];

% Option 2 : For simulated data
valsPerVoxel=WMat{1,1}{1,1};
T2range2plot=T2axis(find(valsPerVoxel(:,1)));

% Plot
figure(5);
for iT2=1:length(T2range2plot)%=50;
	[~ , indT2]=min(abs(T2axis-repmat(T2range2plot(iT2),1,nT2)));
	valsPerVoxel=WMat{1,1}{1,1};
	valsPerVoxel=valsPerVoxel(indT2,:);
	T2Map=im.BW.*(1-im.masked);
	T2Map(indVoxPerSeg{1,:})=valsPerVoxel;
	T2Map2plot=T2Map;
	T2Map2plot(T2Map2plot==0)=NaN;
	subplot(ceil(length(T2range2plot)/2),2,iT2);
	imagesc(T2Map);
	%contourf(X,Y,T2Map2plot);
	colormap hot
	colorbar
	title(sprintf('White Matter T2 = %.f',T2range2plot(iT2)));
	% [X,Y] = meshgrid(1:size((T2Map),2),1:size((T2Map),1));
end


% figure(6);
% imshow(1-im.masked);
% figure(5);
% for iSeg=1:nSegments
%     segMat=cell2mat(WMat{iSeg,1});
%     stem(T2axis,(sum(segMat,2)./sum(sum(segMat,2))),'linestyle','-');grid on;
% 	set(gca,'xscal','log')
% 	hold on;
% end
% axis([T2axis(1) 1e3 0 1.1]);


% Next on the list....
%----------------------------------------------------------------------
% 1) Try different optimization schemes - tikhonov + spatial constraint L1 norm
% 2) Match T2 and T2 distribution to each tube
% 3) DB to try : Jannete fat&water scan + muscle \ Brain scans (manually segmentaiton)
% 4) Mix signals

% 5) Ask Noam about the two options for spatial regulariztion (in one or to iteraitons?)




% Just for testing May 14 [NO]
%sig=0.9*curves(1001,:)'+0.1*curves(1001,:)';
% 		sig=0.8*curves(find(weights(:,3)==1),:)'+0.2*curves(find(weights(:,30)==1),:)';
% 		sigs	= repmat(sig,[1,Nx,Ny]);
% 		sigs	= permute(sigs, [2,3,1] );