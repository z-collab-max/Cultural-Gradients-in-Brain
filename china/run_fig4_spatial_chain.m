function run_fig4_spatial_chain(result_file,spatial_file,out)
% Spatial inputs use the same ROI rows as the China model results.
D=load(result_file);G=load(spatial_file);
if isfield(D,'roi_names'),model_names=D.roi_names;else,model_names=G.model_roi_names;end
assert(isequal(string(model_names(:)),string(G.roi_names(:))),'Model and spatial ROI rows differ');
o=D.behavior_order;tasks=cell(size(o.task_names));tasks(o.best_order)=o.task_names;
effects=D.mediation.task_indirect_z;if ~all(isfinite(effects),'all'),effects=D.mediation.task_indirect;end
for k=1:3
    branches={'behavior','PC1','PC2'};branch=branches{k};
    if k==1,M=D.mediation.behavior_effect;else,M=effects(:,:,k-1);end
    [fit,rho]=best_task_order(M,tasks);order=fit.best_order;F.roi_names=model_names;
    F.task_names=tasks;F.all_orders=fit.orders;F.all_order_rho=rho;F.selected_orders=order;F.selected_order_rho=fit.rsa_map;F.scatter_rsa=fit.rsa_map;
    F.rsa_map=fit.rsa_map;F.surface_slope=fit.slope;F.slope_orders=fit.orders;[~,ranks]=sort(fit.orders,2);x=ranks-mean(ranks,2);F.all_order_slope=M*x'./sum(x.^2,2)';
    if k==1,F.scatter_slope=fit.slope;
    else
        mask=tril(true(size(M,2)),-1);v=fit.model_rdm(mask);v=v-mean(v);F.scatter_slope=(fit.roi_rdm_vectors-mean(fit.roi_rdm_vectors,1))'*v/(v'*v);
    end
    F.ordered_task_effect=M(:,order);F.ordered_task_names=tasks(order);
    if isfield(G,'example_names'),F.example_names=G.example_names(k,:);
    else,[~,hi]=max(fit.slope);[~,lo]=min(fit.slope);F.example_names=G.roi_names([hi lo]);end
    assemble_fig4_spatial(branch,F,G,fullfile(out,branch));
end
end
