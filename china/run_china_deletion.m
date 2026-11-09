function run_china_deletion(input_dir,run_dir,out,progressive,roi_selection)
% Full-sample cultural residuals, fold-specific covariance and artifact removal.
if nargin<4,progressive=true;end
assert(~isfolder(out),'Choose a new output directory');mkdir(out);
[T,Y,center,scores,orders,roi_names]=load_china_run(input_dir,run_dir);
if nargin>=5
    assert(numel(unique(roi_selection))==numel(roi_selection));
    Y=Y(:,roi_selection);roi_names=roi_names(roi_selection);
end
[E,H]=directions(T,Y,[]);
full_clean=center_artifact_endpoint(Y,T,center,scores,orders,H);
full=fit_china_mixed(T,full_clean);full_beta=full.beta;
units=unique(T.province_birth);counts=arrayfun(@(u)sum(T.province_birth==u),units);
[~,order]=sort(counts,'descend');masks={};labels=strings(0,1);kind=strings(0,1);
for k=1:numel(units)
    masks{end+1}=T.province_birth~=units(k);labels(end+1)=string(units(k));kind(end+1)="LOO";
end
if progressive
    for k=1:numel(units)-2
        masks{end+1}=~ismember(T.province_birth,units(order(1:k)));
        labels(end+1)=string(k);kind(end+1)="progressive";
    end
end
fold_beta=nan(size(Y,2),2,numel(masks));status=repmat("rank_deficient",numel(masks),1);
for k=1:numel(masks)
    keep=masks{k};D=T(keep,:);nuisance=[ones(sum(keep),1),D{:,china_covariates}];
    if rank([nuisance,D.culture_pc1,D.culture_pc2])<size(nuisance,2)+2,continue;end
    [retained,~,g]=unique(center(keep));
    if numel(retained)<2,continue;end
    suborders=zeros(numel(retained),size(orders,2));
    for b=1:size(orders,2)
        v=orders(:,b);[ok,indices]=ismember(v,retained);suborders(:,b)=indices(ok);
    end
    [~,H]=directions(D,Y(keep,:),E(keep,:,:));
    clean=center_artifact_endpoint(Y(keep,:),D,g,scores(retained,:),suborders,H);
    fitted=fit_china_mixed(D,clean);fold_beta(:,:,k)=fitted.beta;status(k)="ok";
    rows=find(keep);beta=fitted.beta;save(fullfile(out,sprintf('fold_%03d.mat',k)),'rows','beta','-v7');
end
correlation=nan(numel(masks),2);
for k=find(status=="ok")'
    for pc=1:2,correlation(k,pc)=corr(full_beta(:,pc),fold_beta(:,pc,k));end
end
labels=labels(:);kind=kind(:);
save(fullfile(out,'deletion_source.mat'),'full_beta','fold_beta','correlation','labels','kind','status','roi_names','-v7');
render_china_deletion(fullfile(out,'deletion_source.mat'),out);
end

function [E,H]=directions(T,Y,fixed)
N=height(T);R=size(Y,2);E=zeros(N,R,2);H=E;
for pc=1:2
    name=sprintf('culture_pc%d',pc);
    for r=1:R
        D=T;D.response=zscore(Y(:,r));m=fitlme(D,china_formula('response',{name}),'FitMethod','REML');
        [Psi,sigma2]=covarianceParameters(m);V=sigma2*eye(N);
        groups=m.Formula.GroupingVariableNames;
        for j=1:numel(groups)
            g=findgroups(D.(string(groups{j})));Z=sparse((1:N)',g,1,N,max(g));
            V=V+Psi{j}*full(Z*Z');
        end
        if isempty(fixed)
            X=designMatrix(m,'Fixed');X(:,strcmp(m.Coefficients.Name,name))=[];
            c=D.(name);E(:,r,pc)=c-X*((X'*(V\X))\(X'*(V\c)));
        else
            E(:,r,pc)=fixed(:,r,pc);
        end
        H(:,r,pc)=V\E(:,r,pc);
    end
end
end
