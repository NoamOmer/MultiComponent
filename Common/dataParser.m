

function parsedData=dataParser(dataType)

switch dataType
	
	case 'simulatedData'
		nSegments	= 1;
		Nx			= 10;
		Ny			= Nx;
		type		= 'exprmt';
		disp('You are analyzing simulated data...')
		sim = Simulator;
		sim = sim.getParams(type);			% get parameters for synthesizing experimental EMC
		sim = sim.loadSimData(Dictionary);  % load acquisition operator (E_full)
		sim = sim.simulate(type);			% generate a noisy experimental operator (E_noisy)
		sim = sim.genrtWeights();			% generate vector of weights for the experimental EMC
		sim = sim.genrtExpEMC();			% generate experimental emc
		
		% Create a simulated segment (based on the simulated signal)
		sigs					= repmat(sim.experimental.e,[1,Nx,Ny]);
		sigs					= permute(sigs, [2,3,1] );
		CC.PixelIdxList{:,:}	= find(meshgrid(1:Nx,1:Ny));
		nETL					= size(sigs,3);
		firstTE					= sim.DBprop.firstTE;
		im.masked				= sigs(1,1,1)*ones(Nx,Ny); % Set a simulted image using the values obtained for the first echo
		im.BW					= im.masked;				% In a simulated image the mask and the binary image are the same
		
		% Verify that a correct analysis is performed
		if nETL~=length(TEarr)
			error('The number of echoes in dictionary (TEarr) and simulated data (nETL) missmach !');
		end
		
		if firstTE~=TEarr(1)*1e3 % Convert [uS] to [sec]
			error('The TE in dictionary and TE of the simulated data mismatch !');
		end
		
		
	case 'realData'
		
		% Load MSE data
		%----------------------------------------------------------------------
		analysisType=menu('Type of data to analyze','Dvir''s Phantom','Brain','Fat&Water','MC3','MC3-SingelVoxel');
		
		cd ../../../
		switch  analysisType
			case 1 % Dvir's Phantom
				cd Scans\27_12_17-Dvir\MSE
				analysisType='DvirPhantom';
				vendor			= 'Siemens';
			case 2 % Brain
				cd Scans\Siemens\SEMC149
				analysisType='Brain';
				vendor			= 'Siemens';
			case 3 % Fat&Water
				% 				cd Scans\Siemens\SEMC149
				% 				analysisType='Brain';
				% 				vendor			= 'Siemens';
			case 4
				cd Scans\Bruker\3Comp\2018_05_22\10
				analysisType='MC3';
				vendor			= 'Bruker';
			case 5
				cd Scans\Bruker\3Comp\18_05_22\23
				% Scans\Bruker\3Comp\2018_05_22\24
				analysisType='MC3-SingelVoxel';
				vendor			= 'Bruker';
		end
		SEMCDataDir     = cd;
		
		% 1) Parse SEMCE data
		switch vendor
			case 'Siemens'
				[dcm_info,dcm_info_parsed,im_SEMC_mSl,slices_info] = Read_mSl_Data(SEMCDataDir);
				firstTE  = dcm_info_parsed.alTE(2)*1e-6; % Convert [uS] to [sec]
				nETL	 = size(im_SEMC_mSl,4);
			case 'Bruker'
				st = parse_Bruker_fid_PV5or6([SEMCDataDir,'\'],'fid');
				plotReconImage(st.rdata); % st.rdata -> [nSlices,nETL,Nx,Ny]
				examineSigs(st.rdata,[17,17],1:20); 
				firstTE=st.TE;
				nETL=st.ETL;
				im_SEMC_mSl=permute(st.rdata,[3,4,1,2]);
			case 'GE'
				[dcm_info,dcm_info_parsed,im_SEMC_mSl,slices_info] = Read_mSl_Data_GE(SEMCDataDir);
		end
		
		
		if strcmp(analysisType,'MC3-SingelVoxel')
			% Set the center pixel to be the pixel of interest
			%sigs					= zeros(
			pixOfInterstInd			= [size(im_SEMC_mSl,1),size(im_SEMC_mSl,2)]./2+1;
			nSegments				= 1;
			sigs					= transpose(squeeze(im_SEMC_mSl(pixOfInterstInd(1),pixOfInterstInd(2),:,:))); % [nSlice x nEchos]
			Nx						= 1;
			Ny						= Nx;
			nSlice					= size(sigs,2);
			normSigs				= zeros(Nx,Ny,nSlice,nETL);
			normSigs(1,1,:,:)		= transpose(sigs./repmat(sigs(1,:),nETL,1));
			tmpSigs(1,1,:,:)= sigs';
			sigs					= tmpSigs;
			CC.PixelIdxList{:,:}	= find(meshgrid(1:Nx,1:Ny));
			[indRow,indCol]			= ind2sub(size(squeeze(sigs(:,:,1,1))),CC.PixelIdxList{:,:});
			for iSlice=1:nSlice
				for iETL=1:nETL
					ind(iSlice,iETL)	= sub2ind(size(sigs),indRow,indCol,iSlice,iETL); % [nSlice x nETL]
				end
			end
			
			indVoxPerSeg			= {ind}; % Since only a single voxel is used...
			
		end
		
		nSigsPerSeg				= numel(indVoxPerSeg{:});
end

d.firstTE		= firstTE;
d.nSegments		= nSegments;
d.nSigsPerSeg	= nSigsPerSeg;
d.indVoxPerSeg	= indVoxPerSeg;
d.normSigs		= normSigs;
d.sigs			= sigs;
d.nRows			= Nx;
d.nCols			= Ny;
d.nSlice		= nSlice;
d.nETL			= nETL;

parsedData=d;


% 		% 2) Normalize the measured signals
% 		%	 Henceforth the normalize signals will be analyzed
% 		normSigs	= zeros(size(sigs));
%
% 		switch numel(size(sigs))
% 			case 2 % The single voxel case where sigs -> [nSlice x nEchos]
% 				normSigs = sigs./repmat(sigs(:,1),1,nETL);
% 			case 3 % 1D voxels data case where sigs -> [nRow x nSlice x nEchos]
% 				normSigs	= normalize_SEMC(sigMat);
% 			case 4 % 2D voxels data case where sigs -> [nRow xnCol x nSlice x nEchos]
% 				for iSlice=1:size(sigMat,3)
% 					normSigs(:,:,iSlice,:)= normalize_SEMC(squeeze(sigMat(:,:,iSlice,:)));
% 				end
% 			case 5
% 				normSigs	= normalize_5D_EMC(sigMat);
% 		end
%
% 		sigs		= normSigs;
%
%
% 		% 3) Get the desired ROI
% 		numSlice2plot	= floor(size(sigs,3)/2);	% Select slice for segmentation - currently the middle slice
% 		numTE2plot		= 1;
% 		% In case there is more than one slice
% 		if size(sigMat,4)>1
% 			[im,CC,centers] = getMask(sigs(:,:,numSlice2plot,numTE2plot),analysisType);	% Get mask and segmentation results
% 			CC.PixelIdxList{1,1}(1:20,1)=CC.PixelIdxList{1};
% 			CC.PixelIdxList{1,1}=reshape(CC.PixelIdxList{1,1},5,4);
% 		else % In case there is only one slice
% 			[im,CC,centers] = getMask(sigs(:,:,numSlice2plot),analysisType);	% Get mask and segmentation results
% 		end

% 		nSegments		= size(centers,1);									% Determine the number of segments in the ROI
% 		WMat			= cell(nSegments,1);
% 		candidates		= nan(nCurves,2);
end




