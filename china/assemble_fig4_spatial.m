function assemble_fig4_spatial(branch,F,G,out)
% Spin indices are ROI x rotation, MATLAB one-based. Task null columns match orders.
assert(ismember(branch,{'behavior','PC1','PC2'}));names=string(F.roi_names(:));R=numel(names);
assert(numel(unique(names))==R && isequal(names,string(G.roi_names(:))),'ROI order differs');
labels=string(G.gradient_names(:));expected=["SA_axis";"developmental_expansion";"evolutionary_expansion";"intersubject_variability";"cognitive_PC1";"eigenmode2"];
[ok,col]=ismember(expected,labels);assert(all(ok)&&numel(labels)==6);X=zscore(G.gradients(:,col));
assert(size(X,1)==R && rank(X)==6);spins=double(G.spin_indices);assert(size(spins,1)==R && size(spins,2)>1 && all(spins==fix(spins)&spins>=1&spins<=R,'all'));
assert(all(isfinite([F.all_order_rho,F.selected_order_rho,F.scatter_rsa(:),F.scatter_slope(:),F.surface_slope(:),F.ordered_task_effect]),'all'));
tasks=string(F.task_names(:));T=numel(tasks);orders=double(F.all_orders);assert(T==7 && isequal(sortrows(orders),sortrows(perms(1:T))));
assert(isequal(size(F.all_order_rho),[R size(orders,1)])&&size(F.selected_order_rho,2)==size(F.selected_orders,1));
keep=true(size(orders,1),1);
if strcmp(branch,'behavior')
    groups=["GoNogo","NBack","SDMT","UG","Passball","RO","CVLT"];[ok,j]=ismember(tasks,groups);assert(all(ok));category=[1 1 1 2 2 3 3];cats=category(j);
    keep=~all(cats(orders)==cats(F.selected_orders(1,:)),2);assert(sum(~keep)==24);
end
A.rsa_map=F.rsa_map(:);assert(numel(A.rsa_map)==R && all(isfinite(A.rsa_map)));y=zscore(A.rsa_map);
% Preserve the projection and column order of the original behavior calculation.
if strcmp(branch,'behavior'),sequence=[6 1 2 3 4 5];else,sequence=1:6;end
design=X(:,sequence);project=(design'*design)\design';[~,inverse]=sort(sequence);
taskz=zscore(F.all_order_rho(:,keep));spinz=zscore(y(spins));
beta=project*y;taskbeta=project*taskz;spinbeta=project*spinz;
A.partial_beta=beta(inverse);A.task_order_beta=taskbeta(inverse,:);A.spin_beta=spinbeta(inverse,:);
A.spin_p=mean(abs(A.spin_beta)>=abs(A.partial_beta),2);
A.task_order_p=mean(abs(A.task_order_beta)>=abs(A.partial_beta),2);
target=6;if strcmp(branch,'PC2'),target=1;end
A.pearson_r=X(:,target)'*y/(R-1);A.pearson_task_order_r=X(:,target)'*taskz/(R-1);A.pearson_spin_r=X(:,target)'*spinz/(R-1);
A.pearson_spin_p=mean(abs(A.pearson_spin_r)>=abs(A.pearson_r));
A.pearson_task_order_p=mean(abs(A.pearson_task_order_r)>=abs(A.pearson_r));
A.roi_names=cellstr(names);A.gradient_names=cellstr(expected);A.gradients=G.gradients(:,col);A.task_order_indices_1based=find(keep);A.selected_orders_1based=F.selected_orders;A.all_orders_1based=orders(keep,:);A.task_names=cellstr(tasks);
A.spin_indices_1based=spins;surface_rho=A.rsa_map;
if isfield(F,'surface_rho'),surface_rho=F.surface_rho(:);end
A.surface_z=(surface_rho-mean(F.all_order_rho,2))./std(F.all_order_rho,0,2);
E.scatter_rsa_rho=F.scatter_rsa(:);E.scatter_slope=F.scatter_slope(:);E.surface_slope=F.surface_slope(:);E.roi_names=A.roi_names;
if isfield(F,'slope_color_limits'),E.color_limits=F.slope_color_limits;else,E.color_limits=max(abs(E.surface_slope))*[-1 1];end
E.comparison_r=corr(E.scatter_rsa_rho,E.scatter_slope);E.comparison_spin_null_r=zscore(E.scatter_rsa_rho)'*zscore(E.scatter_slope(spins))/(R-1);
E.spin_exceedance_count=sum(abs(E.comparison_spin_null_r)>=abs(E.comparison_r));E.n_spins=size(spins,2);E.comparison_spin_p=E.spin_exceedance_count/E.n_spins;
E.all_order_slope=F.all_order_slope;E.all_slope_orders_1based=F.slope_orders;E.ordered_task_effect=F.ordered_task_effect;E.ordered_task_names=cellstr(string(F.ordered_task_names));E.spin_indices_1based=spins;
[ok,ix]=ismember(string(F.example_names(:)),names);assert(all(ok));E.example_names=cellstr(names(ix));E.example_task_effect=F.ordered_task_effect(ix,:);
ranks=(1:T)';fit=[ones(T,1),ranks]\E.example_task_effect';E.example_slope=fit(2,:);E.example_r=corr(ranks,E.example_task_effect');
if ~isfolder(out),mkdir(out);end
save(fullfile(out,'fig4_source.mat'),'-struct','A','-v7');save(fullfile(out,'ed12_source.mat'),'-struct','E','-v7');
panels={};panels{1}=struct('folder','Figure4','name',[branch '_gradients'],'kind','violin','beta',A.partial_beta,'p',max(A.spin_p,A.task_order_p),'null',A.task_order_beta','additional_null',A.spin_beta','additional_observed',A.partial_beta,'labels',{A.gradient_names});
panels{2}=struct('folder','Figure4','name',[branch '_alignment'],'kind','scatter','x',A.gradients(:,target),'y',A.rsa_map,'lh',G.left_hemisphere(:),'xlabel',char(expected(target)),'ylabel','RSA Spearman rho','stat_null',A.spin_beta(target,:),'stat_observed',A.partial_beta(target),'stats',sprintf('partial \\beta = %.3f, p_{spin} = %.3f',A.partial_beta(target),A.spin_p(target)));
panels{3}=struct('folder','Figure4','name',[branch '_task_null'],'kind','histogram','null',A.task_order_beta(target,:),'observed',A.partial_beta(target),'p',A.task_order_p(target),'xlabel','Partial coefficient');
for j=1:numel(ix)
    colors=[.70 .49 .45;.44 .51 .62];panels{end+1}=struct('folder','ExtendedDataFigure12','name',sprintf('%s_example%d',branch,j),'kind','task_scatter','x',ranks,'y',E.example_task_effect(j,:)','labels',{E.ordered_task_names},'color',colors(1+mod(j-1,2),:),'stats',sprintf('slope = %.4f, r = %.3f',E.example_slope(j),E.example_r(j)));
end
panels{end+1}=struct('folder','ExtendedDataFigure12','name',[branch '_rsa_slope'],'kind','identity_scatter','x',zscore(E.scatter_rsa_rho),'y',zscore(E.scatter_slope),'xlabel','RSA Spearman rho (z)','ylabel','Slope (z)','fit_line',true,'stats',sprintf('r = %.3f, spin exceedances = %d/%d',E.comparison_r,E.spin_exceedance_count,E.n_spins));
save(fullfile(out,'panels.mat'),'panels','-v7');
end
