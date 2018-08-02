
function [weights] = generateWeights(p)

% Get parameters


minHeight   = p.minHeight;      % [#]
maxHeight   = p.maxHeight ;     % [#]
dH          = p.dH;             % [#]
nT2         = p.nT2;            % [#]
isSave		= 1;			    % save flag - '0' by default (set to '1' to save weights to a selected directory)	
nMaxComp    = p.nMaxComp ;      % Number of components within dictionary
counter     = 0;                % Weights counter
T2axis      = p.T2axis;
nHeights    = length(minHeight:dH:maxHeight);
heights     = linspace(minHeight,maxHeight,nHeights);
nIter 		= ceil((factorial(nT2)/(factorial(nMaxComp)*factorial(nT2-nMaxComp)))* ...
			  ((nHeights)*(nHeights+1)/2));
weights     = zeros(length(T2axis),nIter);


h_wb = waitbar(0,'Generating weights DB ...');
tmp2=1;
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
						
						counter	= counter+1;
						wv		= zeros(length(T2axis),1);
						wv(i)	= heights(hi);
						wv(j)	= heights(hj);
						wv(k)	= heights_hk;
						
						weights(:,counter)=wv;
						waitbar(counter / nIter,h_wb,...
						sprintf('Generating weights DB -  %d / %d',counter,nIter));

					end
					
				end
				
			end
			
		end
		
	end
	
end

% Calculate the double component case for the last 2 positions 
i=(nT2-1);
for hi= nHeights:-1:1
	
	wv					= zeros(length(T2axis),1);
	wv(i)				= heights(hi);
	wv(i+1)				= maxHeight-heights(hi);
	counter				= counter+1;
	weights(:,counter)	= wv;
	
end

% Make sure uniqueness of all weights
weights = unique(weights','rows');
weights = weights';
if sum(weights(:,1))== 0
	weights(:,1)=[];
end

close(h_wb);

% Save the weights permutations matrix

if isSave
    
    % Option 1 - save to the 'WeightsPermutations' folder (under results) 
    orgCD=cd;
    compName=getComputerName();
    setCD(compName);
    cd ..
    cd results\WeightsPermutations
    save(['w_T2_',num2str(nT2),'_H_',num2str(nHeights)],'weights','p');
    cd(orgCD);
    
    % Option 2 - save to a selected folder 
    %folderName=uigetdir(pwd,'Select a directory to save the weights matrix');
    %save([folderName,'\','w_T2_',num2str(nT2),'_H_',num2str(nHeights)],'weights','p');
end

% Display notification to command window
display('Weights permutations matrix was generated !');




