function result=run_sensitivity(T,cortex,additional_sets)
% additional_sets is a cell array of covariate-name cell arrays.
T=prepare_china_table(T);result=cell(size(additional_sets));
for j=1:numel(additional_sets)
 names=additional_sets{j};D=T;
 for k=1:numel(names),D.(names{k})=zscore(D.(names{k}));end
 result{j}=fit_china_mixed(D,cortex,names);
end
end
