
function [im,CC,centers]=getMask(img,scanType)


% Configuration
switch scanType
    
    case 'DvirPhantom'
        nSegs					= 12; %6
        tubeDiameter			= 25; % [pixels] in the future change to cm
        seErodeSize				= 2;
        seOpenSize				= 3;
        factor4excludingSegment = .2;
        
    case 'SPIOPhantom'
        nSegs					= 6;
        tubeDiameter			= 28;
        seErodeSize				= 2;
        seOpenSize				= 3;
        factor4excludingSegment = .2;
        
    case 'SPIOPhantom_2nd'
        nSegs					= 7;
        tubeDiameter			= 8;
        winsiz                  = floor(tubeDiameter/2);
        win                     = -winsiz:winsiz;
        seErodeSize				= 2;
        seOpenSize				= 3;
        factor4excludingSegment = .2;
        % Select the centers of each tube manually and craete the mask accordinglly
        %       [200[ug/ml]    150[ug/ml]   100[ug/ml]     50[ug/ml]   200[ug/ml]- from 1st batch   MnCl    Water ]
        centers=[37,71 ;        41,40 ;      91,32 ;        106,69 ;    72,104 ;                   73,74 ; 65,51];
        
        for iTube =1:nSegs
            img(centers(iTube,2)+win,centers(iTube,1)+win)=100.*ones(size(win,2),size(win,2));
        end
        img(img~=100)=0;
        img(img==100)=1;
        CC=[];
        im=img;
        return;
        
    case 'SPIOPhantom_3nd'
        nSegs					= 8;
        tubeDiameter			= 10;
        winsiz                  = floor(tubeDiameter/2);
        win                     = -winsiz:winsiz;

        % Select the centers of each tube manually and craete the mask accordinglly
        %       [75[ug/ml]  50[ug/ml]   25[ug/ml]   10[ug/ml]   5[ug/ml]    2.5[ug/ml]  MnCl      Water ]
        centers=[96,103 ;   139,54 ;    162,98 ;    143,161 ;   71,167 ;    28,126 ;    32,77   ; 81,34];
        
        maskArbValue=1e9;
        for iTube =1:nSegs
            img(centers(iTube,2)+win,centers(iTube,1)+win)=maskArbValue.*ones(size(win,2),size(win,2));
        end
       
        img(img~=maskArbValue)  =   0;
        img(img==maskArbValue)  =   1;
        CC                      =   [];
        im=img;
        return;
        
    case 'Brain' %Only WM at this stage
        seCloseSize				= 5;
        im.org          = img;
        im.median       = medfilt2(im.org,[2 3]);
        tmp             = im.median;
        tmp(tmp==0)     = 999;
        tmp(tmp~=999)   = 0;
        se2             = strel('disk',seCloseSize);
        se3             = strel('disk',12);
        
        im.BWclose      = imclose(tmp,se2);
        tmp2            = edge(im.BWclose);
        CC              = bwconncomp(1-im.BWclose);
        mask            = regiongrowing(1-im.BWclose,53,69,0.9);
        im.masked=im.org; im.masked(~mask)=0; % imshow(im.masked)
        im.masked=1-im.masked;
        
        [counts,binLocations]	= imhist(im.masked);					% Calculate histogram for the stretched image
        nPixels					= sum(counts);							% Count the number of pixels within the image
        counts(end)				= 0;									% Exclude the backgorund (space+phantomo boundaries) from histogram
        counts(1)               = 0;
        counts(end)             = 0;
        hitogramThr				= otsuthresh(counts);					% Find threhsold based on stretched histogram (Otsu's)
        im.BW					= imbinarize(im.masked,0.55);	% Convert to a binary image based on threshold
        im.BW(~mask)=0;
        maskErode     = imerode(mask,se3);
        im.BW					= im.BW.*maskErode;
        CC                      = bwconncomp(im.BW);
        numPixels               = cellfun(@numel,CC.PixelIdxList);
        [~,indNumPixelsSorted]  = sort(numPixels,'descend');
        
        CC.PixelIdxList={CC.PixelIdxList{1,(indNumPixelsSorted(1))}};
        CC.NumObjects=1;
        %centers = regionprops(CC,'Centroid');
        
        h_figure = figure();
        screensize = get( 0, 'Screensize' );
        screensize(1)=10;screensize(2)=60;
        h_figure.Position = .8.*screensize;
        subplot(3,3,1);imshow(im.org);							title('Original');
        subplot(3,3,4);imshow(mask);							title('Mask');
        subplot(3,3,[2 3]);plot(binLocations,counts./nPixels);grid; title(sprintf('Histogram Threshold %.2f',hitogramThr));
        subplot(3,3,7);imshow(im.masked);							title('Image * Mask ');
        subplot(3,3,[5 6 8 9]);imshow(im.BW);				title('Final Mask');
        
        choice = questdlg('Are you pleased with segentation results ?','Select an answer');
        close(h_figure);
        centers = [];
        switch choice
            case 'Cancel'
                im = [];
                CC = [];
                disp('Segmentation was canceled');
                return;
            case 'No'
                disp('Retry to perform segmentation with different paramenters');
                im = [];
                CC = [];
                return;
                
            case 'Yes'
                % Save the centers of each component
                centers = regionprops(CC,'Centroid');
                return;
        end
        
        return;
        % Select ROIs manually
        %curRoi= roipoly(im.org);
        %[vall indd]=max(curRoi);
        %roi{iSeg} =curRoi;
        %CC.PixelIdxList{1,iSeg}=find(roi{1,iSeg});
        
		
		
		case 'MC3-SingelVoxel'
        nSegs					= 1;
        tubeDiameter			= 1;
        winsiz                  = floor(tubeDiameter/2);
        win                     = -winsiz:winsiz;

        % All 3 tubes are within a single voxel located at the center of the image 
        centers=[17,17];
        
        maskArbValue=1e9;
        for iTube =1:nSegs
            img(centers(iTube,2)+win,centers(iTube,1)+win)=maskArbValue.*ones(size(win,2),size(win,2));
        end
       
        img(img~=maskArbValue)  =   0;
        img(img==maskArbValue)  =   1;
        CC.PixelIdxList{1,1} =  sub2ind(size(img),centers(1),centers(2)) 
        im=img;
        return;
		
        
end

% Allocations
se1							= strel('disk',seErodeSize);
se2							= strel('disk',seOpenSize);
estimation4nPixsPerSeg		= .25*pi*tubeDiameter^2;				% Assuming circular tubes

% Process
im.org						= img;									% Extract the middle slice
im.median					= medfilt2(im.org,[2 3]);				% Clean salt&pepper noise (assuming smoothness)
im.spacePix					= find(im.median>=1);					% Store the space pixels indexes
im.phantomBounds			= find(im.median==0);					% Store the phantom-air pixels indexes
im.median(im.phantomBounds) = 1;									% Change phantom boundaries to background
im.stretch					= mat2gray(im.median,[0 1]);			% Stretch the image intensities to be in the range of [0 1]
[counts,binLocations]		= imhist(im.stretch);					% Calculate histogram for the stretched image
nPixels						= sum(counts);							% Count the number of pixels within the image
counts(end)					= 0;									% Exclude the backgorund (space+phantomo boundaries) from histogram
hitogramThr					= otsuthresh(counts);					% Find threhsold based on stretched histogram (Otsu's)
im.BW						= imbinarize(im.stretch,hitogramThr);	% Convert to a binary image based on threshold
im.BW						= 1-im.BW;								% Replace background and object
im.BWerode					= imerode(im.BW,se1);					% Erode image to remove small objectes (see se1 size)
im.BWopen					= imopen(im.BWerode,se2);				% Open image to smooth tubes edges
CC							= bwconncomp(im.BWopen);				% Find connected components in binary image
numPixels					= cellfun(@numel,CC.PixelIdxList);		% Create a vector of the number of elements per segment
im.final					= im.BWopen;							% Store the final processed image under 'final' field

[compSortedPerSize,compSortedPerSizeInd]=sort(numPixels,'descend');
indComp2exclude =	find (	compSortedPerSize	>	(1+factor4excludingSegment)*estimation4nPixsPerSeg |...
    compSortedPerSize	<	(1-factor4excludingSegment)*estimation4nPixsPerSeg		);

%Exclude components which are not tubes ( thier size deviates from thier expected size)
for i_comp=1:length(indComp2exclude)
    im.final(CC.PixelIdxList{indComp2exclude(i_comp)}) = 0;
    CC.PixelIdxList(indComp2exclude(i_comp))=[];
end

% Return in case not all tubes were found \ too many components were found
CC.NumObjects=CC.NumObjects-length(indComp2exclude);
if nSegs~=(size(CC.PixelIdxList,2)-indComp2exclude)
    disp('The nunber of segments don''t match to the number of phantom tubes !')
    return;
end


% Plot processed images per stage
h_figure = figure();
screensize = get( 0, 'Screensize' );
screensize(1)=10;screensize(2)=60;
h_figure.Position = .8.*screensize;
subplot(3,3,1);imshow(im.org);							title('Original');
subplot(3,3,4);imshow(im.stretch);						title('Median & Stretched');
subplot(3,3,[2 3]);plot(binLocations,counts./nPixels);grid; title(sprintf('Histogram Threshold %.2f',hitogramThr));
subplot(3,3,7);imshow(im.BW);							title('Binary');
subplot(3,3,[5 6 8 9]);imshow(im.BWopen);				title('Final Mask');

choice = questdlg('Are you pleased with segentation results ?','Select an answer');
close(h_figure);
centers = [];
switch choice
    case 'Cancel'
        im = [];
        CC = [];
        disp('Segmentation was canceled');
        return;
    case 'No'
        disp('Retry to perform segmentation with different paramenters');
        im = [];
        CC = [];
        return;
        
    case 'Yes'
        % Save the centers of each component
        centers = regionprops(CC,'Centroid');
        return;
end

end