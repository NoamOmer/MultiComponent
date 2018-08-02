function []= plotReconImage(reconImage,echos2plot)
% This function plots the reconstructed images for a given set of data

[nSlice,nETL,Nx,Ny]=size(reconImage);

if nargin<2
	echos2plot=[1,4,nETL];
end
for iETL=echos2plot
	figure(1);
	
	for iPlot=1:nSlice
		% Make sure all signals maps are presented with the same axes
		subplot(4,ceil(nSlice/4),iPlot);
		imagesc(squeeze(reconImage(iPlot,iETL,:,:)));
		xlabel(['Slice: ',num2str(iPlot)])
		colorbar;
		if (iETL==1 && iPlot==1);		cl = caxis;		end
		caxis(cl);
		grid on;
	end
	mtit(['Reconstructed Images - Echo: ',num2str(iETL)]);
	k = waitforbuttonpress;
end
