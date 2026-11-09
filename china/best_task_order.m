function [result,all_rho]=best_task_order(effects,task_names)
% effects: ROI-by-task coefficients; signed lower-triangle RDMs.
[n_rois,n_tasks]=size(effects);orders=sortrows(perms(1:n_tasks));mask=tril(true(n_tasks),-1);
roi_vectors=zeros(nnz(mask),n_rois);models=zeros(nnz(mask),size(orders,1));
for r=1:n_rois,v=effects(r,:)';D=v-v';roi_vectors(:,r)=D(mask);end
for j=1:size(orders,1)
 ranks=zeros(n_tasks,1);ranks(orders(j,:))=(1:n_tasks)';D=ranks-ranks';models(:,j)=D(mask);
end
roi_rank=zscore(tiedrank(roi_vectors));model_rank=zscore(tiedrank(models));
all_rho=roi_rank'*model_rank/(nnz(mask)-1);assert(all(isfinite(all_rho),'all'));
score=mean(abs(all_rho),1);best=find(score>=max(score)-1e-12,1);
ranks=zeros(n_tasks,1);ranks(orders(best,:))=(1:n_tasks)';
result=struct('best_order',orders(best,:),'task_names',{task_names(orders(best,:))},...
 'score',score(best),'rsa_map',all_rho(:,best),'all_fit_scores',score,...
 'orders',orders,'model_rdm',ranks-ranks','roi_rdm_vectors',roi_vectors,...
 'slope',((ranks-mean(ranks))'*(effects-mean(effects,2))')'/sum((ranks-mean(ranks)).^2));
end
