
clcl;

% pcName=getComputerName;
% setCD(pcName);
%fp='C:\Users\me\Google Drive\Scans\27_12_17-Dvir\MSE';

fp = 'C:\Users\Itzhak\Google Drive\Scans\27_12_17-Dvir\MSE';

data                = dcmReadFrames(fp,'MSE');
[nRows,nCols,nEchos]= size(data);
nVox                = nRows*nCols;
e                   = cell(nVox,1);
rowInd              = zeros(nRows,1);
colInd              = zeros(nCols,1);
voxCounter          =0;

for i_rows=1:nRows
	for i_cols=1:nCols
		voxCounter          =voxCounter+1;
		e{voxCounter,:}     = squeeze(data(i_rows,i_cols,:));
        rowInd(voxCounter)  = i_rows;
        colInd(voxCounter)  = i_cols;
	end
end

dataTable = table(rowInd,colInd,e);
ROI= [40 50 100 110]; % [rows:buttom(1)->top(2) cols:buttom(3)->top(4)
%figure();imagesc(data(:,:,30));

% Extract signal form ROI 
ROInd = dataTable.rowInd >= ROI(1) &  dataTable.rowInd < ROI(2) & ...
       dataTable.colInd >= ROI(3) &  dataTable.colInd < ROI(4) ;
  
e     = reshape(cell2mat(dataTable.e(ROInd,1)),[],sum(ROInd));


figure();plot(e,'.-');grid;
   
