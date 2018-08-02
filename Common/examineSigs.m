function examineSigs(sigData,vox2plot,slices)
% This function plots the acquired signal for a given voxel 

[nSlice,nETL,Nx,Ny]=size(sigData);

if nargin<2
	% By default plot the pixel at the center of the sample
	vox2plot=[floor(Nx/2)+1,floor(Ny/2)+1];
	slices=floor(nSlice/2)+1;
end

figure(1);
for iSlice=1:nSlice
	for iNx = vox2plot(1)
		for iNy = vox2plot(2)
			sig2Plot=squeeze(sigData(iSlice,:,iNx,iNy))';
			subplot(2,1,1);	plot(sig2Plot./repmat(sig2Plot(1,:)',1,nETL)','.-'); grid on; ylabel('Normalized Signal [A.U]');xlabel('Echo [#]');
			subplot(2,1,2);	plot(sig2Plot,'.-');grid on;ylabel('Signal [A.U]');xlabel('Echo [#]');
			mtit(['Voxel : [',num2str(iNx),'/',num2str(Nx),' ,',num2str(iNy),'/',num2str(Nx),' Slice: ',num2str(iSlice),']']);
			k = waitforbuttonpress;
		end
	end
end
		
end



