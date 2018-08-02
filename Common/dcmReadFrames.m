function data=dcmReadFrames(fp,type)

data=[];

if strcmp(type,'MSE')
	files=dir([fp,'/*.dcm']);
	nDcmFiles=size(files,1);
	
	%Rows: 128
	%Columns: 128
	%data(:,:,i_dcmFile)=dicominfo(fullfile(fp,files(1,i_dcmFile).name));
	
	for i_dcmFile=1:nDcmFiles
		data(:,:,i_dcmFile)=dicomread(fullfile(fp,files(i_dcmFile,1).name));
	end
	
end
	
end