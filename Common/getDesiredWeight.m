
function desiredWeightInd=getDesiredWeight(weights,locs,relFract,T2axis)
% Using eps is required since from some reason some of the weights are not accurate 
% One way of fixing it is to round all permutations before the database creation

A	= find(weights(:,find(T2axis==locs(1)))>=(relFract(1)-eps) & weights(:,find(T2axis==locs(1)))<=(relFract(1)+eps));
B	= find(weights(:,find(T2axis==locs(2)))>=(relFract(2)-eps) & weights(:,find(T2axis==locs(2)))<=(relFract(2)+eps));
C	= find(weights(:,find(T2axis==locs(3)))>=(relFract(3)-eps) & weights(:,find(T2axis==locs(3)))<=(relFract(3)+eps));
AB	= intersect(A,B);

desiredWeightInd	= intersect(AB,C);

if isempty(desiredWeightInd)
    desiredWeightInd=[];	disp('Dictionary does not include the desired weight');
end

end