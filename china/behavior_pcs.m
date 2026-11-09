function [scores,metadata]=behavior_pcs(behavior,directions)
if nargin<2,directions=cellfun(@(x)ones(1,size(x,2)),behavior,'UniformOutput',false);end
scores=cell(size(behavior));metadata=cell(size(behavior));
for t=1:numel(behavior)
 X=zscore(behavior{t}.*directions{t});
 assert(all(isfinite(X),'all'),'Complete behavioral inputs are required for PCA.');
 [coeff,score,~,~,explained]=pca(X);k=find(cumsum(explained)>85,1);
 if isempty(k),k=size(score,2);end
 signs=sign(sum(coeff(:,1:k),1));
 for j=find(abs(sum(coeff(:,1:k),1))<1e-12)
  [~,ix]=max(abs(coeff(:,j)));signs(j)=sign(coeff(ix,j));
 end
 scores{t}=zscore(score(:,1:k).*signs);
 metadata{t}=struct('loadings',coeff(:,1:k).*signs,'explained_percent',explained,'retained',k,'metric_directions',directions{t});
end
end
