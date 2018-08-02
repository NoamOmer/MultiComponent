clcl;
%% Params
nComponent  = 2;	% [#]
minT2       = 1;	% [ms]
maxT2       = 1e3;	% [ms]
dT2         = 1;	% [ms]
nHeights    = 10;	% [#]
minHeight   = 1;	% [#]
maxHeight   = 10;	% [#]
nWidth      = 11;	% [#]
minWidth    = 5;	% [#]
maxWidth    = 15;	% [#]
nT2         = 30;	% [#]
min_T2_dist = 3;	
T2base      = 30;	% [ms]
round_f     = 1;

%% 
T2axisFull  = minT2:dT2:maxT2;  
heights     = linspace(minHeight,maxHeight,nHeights);
widths      = linspace(minWidth,maxWidth,nWidth);
T2axis      = calc_equispaced_T2(minT2,maxT2,T2base,nT2,min_T2_dist,round_f,0);
weights     = zeros(length(T2axisFull),nT2*nHeights*nWidth*sum(1:nComponent));

nIterPerComponent	= length(T2axis)*nHeights*nWidth;
nIter				= nIterPerComponent.^(nComponent);

mu		= zeros(nIter,1);
sigma	= zeros(nIter,1);
amp		= zeros(nIter,1);
i_w		= 0; 
k		= 1;


h_wb = waitbar(0,'Generating weights DB ...');
figure();

for i_comp = 1:nIter
 
for i_T2 = 1:length(T2axis)
for i_ht = 0:(nHeights-1)  
for i_wd = 1:nWidth 

    i_w=i_w+1;
	waitbar(i_w / nIter,h_wb,sprintf('Generating weights DB -  %d / %d',i_w,nIter))
    
    mu(i_w)    = T2axis(i_T2);
    sigma(i_w) = widths(i_wd);
    amp(i_w)   = heights(i_ht+1);
    
    if ~mod(i_w,nIterPerComponent);        k=k+1;    end
    
	weights = amp(i_w)*gaussmf(T2axisFull,[sigma(i_w) mu(i_w)])'; % weights(:,k)+
    %weights(:,i_w)= weights(:,k)+amp*gaussmf(T2axisFull,[sigma mu])';
 
    % Plot 
    plot(weights); 
    axis([0 maxT2 0 10]);grid;
	title(sprintf('\\mu: %.f [mS] \\sigma: %.f [mS] Amp: %.f ',mu(i_w),sigma(i_w),amp(i_w)))
    pause(1e-5); 
	
    
end
end
end
end


close(h_wb) 



