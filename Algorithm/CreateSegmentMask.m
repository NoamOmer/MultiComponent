
function mask_im = createSegmentMask(im,segType)

if nargin<2;	segType = 'WM'; end

mask_im   = ones(size(im));
horizontalLimits = [50:100];
verticalLimits = [50:100];
segmetSize=length(horizontalLimits)*length(verticalLimits);
mask_im(verticalLimits,horizontalLimits)=0;
figure();imshow(mask_im);
end