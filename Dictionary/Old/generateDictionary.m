clcl;

tic

%% Params
minT2       = 1;	% [ms]
maxT2       = 1e3;	% [ms]
dT2         = 1;	% [ms]
minHeight   = 0;	% [#]
maxHeight   = 6;	% [#]
dH			= 1;	% [#]
nT2         = 25;	% [#]
min_T2_dist = 2;	% [ms]
T2base      = 40;	% [ms]
round_f     = 1;    % plot-flag {'0','1'} (set to '1' to plot T2 axis)
nMaxComp    = 3;	% Number of components within dictionary;

%%
% T2axisFull  =	minT2:dT2:maxT2;
nHeights    =	length(minHeight:dH:maxHeight);
heights     =	linspace(minHeight,maxHeight,nHeights);
T2axis      =	calc_equispaced_T2(minT2,maxT2,T2base,nT2,min_T2_dist,round_f,0);
nIter 		=	ceil((factorial(nT2)/(factorial(nMaxComp)*factorial(nT2-nMaxComp)))* ...
    ((nHeights)*(nHeights+1)/2));
weights     =	zeros(length(T2axis),nIter);
stamp		=	zeros(6,nIter);

counter=0;
h_wb = waitbar(0,'Generating weights DB ...');
%a=0;b=0;c=0;d=0;e=0;tmp_j=0;tmp_i=0;
%figure;

tmp2=1;
tmp3=1;
for i = 1:(nT2-(nMaxComp-1))
    tmp=1;
    
    for j = (i+1):(nT2-(nMaxComp-2))
        
        if j>(i+2)
            tmp2=1;
        end
        
        for k = (j+1):(nT2)
            
            for hi=nHeights:-1:2
                
                if (hi == nHeights)
                    if tmp==1
                        counter=counter+1;
                        wv = zeros(length(T2axis),1);
                        wv(i) = heights(hi);
                        weights(:,counter)=wv;
                        tmp=0;
                        continue;
                    else
                        continue;
                    end
                end
                
                for hj=1:(nHeights+1-hi)
                    
                    sumHeights = heights(hi)+ heights(hj);
                    heights_hk = maxHeight-sumHeights;
                    
                    if  ((sumHeights+heights_hk)== maxHeight)  && (hj<=nHeights) && (heights_hk<=maxHeight)
                        
                        if (sumHeights == maxHeight)
                            if tmp2<(nHeights-1)
                                counter=counter+1;
                                %stamp(counter,1:6)=[i,j,k heights(hi), heights(hj), heights_hk];
                                wv = zeros(length(T2axis),1);
                                wv(i) = heights(hi);
                                wv(j) = heights(hj);
                                wv(k) = heights_hk;
                                weights(:,counter)=wv;
                                tmp2=tmp2+1;
                                continue;
                            else
                                continue;
                            end
                        end
                        
                        
                        counter=counter+1;
                        % stamp(counter,1:6)=[i,j,k heights(hi), heights(hj), heights_hk];
                        wv = zeros(length(T2axis),1);
                        wv(i) = heights(hi);
                        wv(j) = heights(hj);
                        wv(k) = heights_hk;
                        
                        weights(:,counter)=wv;
                        
                        %stem(wv);grid on;
                        %axis([0, 50, 0, 10]);
                        %pause(0.0001);
                        waitbar(counter / nIter,h_wb,...
                            sprintf('Generating weights DB -  %d / %d',counter,nIter));
                        
                        % Save the weights matrix from time to time
                        if (mod(counter,1e5)==0)
                            save('weights','weights');
                        end
                        
                    end
                    
                end
            end
            
        end
        
    end
end

% Include the weights of the last 2 locations
i=(nT2-1);
for hi= nHeights:-1:1
    
    wv = zeros(length(T2axis),1);   
    wv(i) = heights(hi);
    wv(i+1) = maxHeight-heights(hi);
    counter=counter+1;
    weights(:,counter)=wv;
end

% Make sure the dictionary uniqueness
weights=unique(weights','rows');
weights=weights';
if sum(weights(:,1))== 0
    weights(:,1)=[];
end

% Save Results
save(['w_T2_',num2str(nT2),'_H_',num2str(nHeights)],'weights');
runTime=toc;





%sumHeights = heights(hi)+ heights(hj);
%if sumHeights==maxHeight
%end

% 						if (j+1==k && heights(hj)==0)
% 							if tmp3<(nHeights-1)
% 								wv = zeros(length(T2axis),1);
% 								wv(i) = heights(hi);
% 								wv(j) = heights(hj);
% 								wv(k) = heights_hk;
% 								tmp3=tmp3+1;
% 								continue;
% 							else
% 								continue;
% 							end
% 						end
%