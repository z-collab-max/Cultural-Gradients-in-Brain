function [shuffled,order]=shuffle_china_culture(T,order)
% Shuffle province-level culture rows jointly across cultural PCs.
groups=unique(T.province_birth,'stable');G=numel(groups);values=zeros(G,2);
for j=1:G
 hit=T.province_birth==groups(j);
 values(j,:)=[T.culture_pc1(find(hit,1)),T.culture_pc2(find(hit,1))];
 assert(all(abs(T.culture_pc1(hit)-values(j,1))<1e-10));
 assert(all(abs(T.culture_pc2(hit)-values(j,2))<1e-10));
end
if nargin<2,order=randperm(G);end
assert(isequal(sort(order),1:G));shuffled=T;
[present,birth]=ismember(T.province_birth,groups);assert(all(present));
[present,destination]=ismember(T.province_destination,groups);assert(all(present));
for pc=1:2
 x=sprintf('culture_pc%d',pc);y=sprintf('destination_pc%d',pc);
 v=values(order,pc);shuffled.(x)=v(birth);shuffled.(y)=v(destination);
end
end
