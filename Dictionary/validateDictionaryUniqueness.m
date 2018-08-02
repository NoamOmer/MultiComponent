
%load('weights.mat')

% Clear all non-weights 
% (pure zeros vector generate as a results of memory allocation)
indZero=zeros(size(weights,2),1);
c=0;
for iw=1:size(weights,2)
	if sum(weights(:,iw))==0
		c=c+1;
		indZero(c)=iw;
	end
end

indZero(indZero==0)=[];
weights(:,indZero)=[];


% Extract a stamp for each weight from the dictionary
stamp=[];
stamp=cell(size(weights,2),1);
for iw=1:size(weights,2)
	curW=weights(:,iw);
	curInd=find(curW);
	curVal=curW(curInd);	
	stamp{iw}=[curInd ; curVal]';
end

i2=0;i4=0;i6=0;stamp2=[];stamp4=[];stamp6=[];
for ii=1:numel(stamp)
	switch length(stamp{ii})
		case 2
			i2=i2+1;
			stamp2(i2,1:2)=stamp{ii};
		case 4
			i4=i4+1;
			stamp4(i4,1:4)=stamp{ii};
		case 6
			i6=i6+1;
			stamp6(i6,1:6)=stamp{ii};
	end
end

notUnique=0;
for ii=1:size(stamp2,1)
	for jj=1:size(stamp2,1)
	if sum(stamp2(ii,:)==stamp2(jj,:))==size(stamp2,2) & (jj~=ii)
		notUnique=ii;
	end
	end
end

for ii=1:size(stamp4,1)
	for jj=1:size(stamp4,1)
	if sum(stamp4(ii,:)==stamp4(jj,:))==size(stamp4,2) & (jj~=ii)
		notUnique=notUnique+1;
	end
	end
end


for ii=1:size(stamp6,1)
	for jj=1:size(stamp6,1)
	if sum(stamp6(ii,:)==stamp6(jj,:))==size(stamp6,2) & (jj~=ii)
		notUnique=ii;
	end
	end
end


x=1;

% % Roger Stafford's solution
% B1 = randfixedsum(nExprmnt,nT2,1,0,1);
% % "Easy method" [Please don't use!]
% A = rand(100000,10);
% B2 = bsxfun(@rdivide,A,sum(A,2));
% % Compare the marginal distributions of one variable
% figure
% subplot(2,1,1),hist(B1(1,:),25)
% subplot(2,1,2),hist(B2(:,1),25)
% 
% 
% 
% options = optimoptions(@lsqlin,'MaxIter',10000);
% w_fitted = lsqlin(E,e,[],[],[],[],[],[],[],options);