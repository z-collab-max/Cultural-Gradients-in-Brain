function run_dimension_branches(input_dir,run_dir,dimension_file,out,inner_draws,seed,roi_selection)
% Fixed full-sample ComBat; predictor-specific artifact endpoints and mixed models.
if nargin<5,inner_draws=10000;end
if nargin<6,seed=42;end
assert(~isfolder(out),'Choose a new output directory');mkdir(out);
[T,Y,center,~,~,roi_names]=load_china_run(input_dir,run_dir);
if nargin>=7,Y=Y(:,roi_selection);roi_names=roi_names(roi_selection);end
U=readtable(dimension_file,'TextType','string');assert(strcmp(U.Properties.VariableNames{1},'province'));
assert(numel(unique(U.province))==height(U));dimension_labels=U.Properties.VariableNames(2:end);V=U{:,2:end};
assert(size(V,2)>=3&&all(isfinite(V),'all')&&all(std(V)>0));
[found,birth]=ismember(string(T.province_birth),string(U.province));assert(all(found));
centers=splitapply(@(x)x(1),string(T.center_province),center);
[found,ci]=ismember(centers,string(U.province));assert(all(found));
stream=RandStream('mt19937ar','Seed',seed);orders=zeros(numel(ci),inner_draws);
for b=1:inner_draws,orders(:,b)=randperm(stream,numel(ci));end
culture=readtable(fullfile(input_dir,'culture.csv'),'TextType','string');[found,ui]=ismember(string(U.province),string(culture.province));assert(all(found));
base=culture{ui,{'pc1','pc2'}};
clean=center_artifact_endpoint(Y,T,center,base(ci,:),orders);fit=fit_china_mixed(T,clean);full_maps=fit.beta';
R=size(Y,2);D=size(V,2);dimension_effects=zeros(R,D);fold_maps=zeros(D,2,R);matching=zeros(D,2);orientation=matching;
for omit=1:D
 C=[V(:,omit),V(:,omit)];Q=replace_scores(T,C,birth);
 clean=center_artifact_endpoint(Y,Q,center,C(ci,:),orders);fit=fit_china_mixed(Q,clean);dimension_effects(:,omit)=fit.beta(:,1);
 [~,scores]=pca(zscore(V(:,setdiff(1:D,omit))));r=corr(base,scores(:,1:2));
 if sum(abs(diag(r)))>=sum(abs(diag(r(:,[2 1])))),idx=[1 2];else,idx=[2 1];end
 sg=[sign(r(1,idx(1))),sign(r(2,idx(2)))];assert(all(sg~=0));C=scores(:,idx).*sg;
 Q=replace_scores(T,C,birth);clean=center_artifact_endpoint(Y,Q,center,C(ci,:),orders);fit=fit_china_mixed(Q,clean);
 fold_maps(omit,:,:)=fit.beta';matching(omit,:)=idx;orientation(omit,:)=sg;
 fprintf('Completed cultural dimension %d/%d\n',omit,D);
end
save(fullfile(out,'dimension_models.mat'),'full_maps','dimension_effects','fold_maps','matching','orientation','dimension_labels','roi_names','inner_draws','seed','-v7');
end
function Q=replace_scores(T,C,birth)
Q=T;
for pc=1:2,Q.(sprintf('culture_pc%d',pc))=zscore(C(birth,pc));end
end
