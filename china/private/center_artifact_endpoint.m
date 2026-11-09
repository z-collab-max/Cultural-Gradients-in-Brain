function [clean,kept]=center_artifact_endpoint(Y,T,center,center_scores,orders,protection)
% Center-level score shuffles, ROI-specific OLS contributions and PCA residuals.
counts=accumarray(center,1);n=numel(center);k=numel(counts)-1;
X=[ones(n,1),T{:,china_covariates}];[Q,~]=qr(X,0);
assert(rank(X)==size(X,2),'Singular artifact nuisance design');
Y=zscore(Y);clean=zeros(n,size(Y,2),2);kept=false(size(orders,2),2);
for pc=1:2
    target=T.(sprintf('culture_pc%d',pc));scores=center_scores(orders,pc);
    scores=reshape(scores,size(orders));
    % Correlation with individual culture, weighted by the number at each center.
    mu=counts'*scores/n;centered=scores-mu;
    scale=sqrt(counts'*(centered.^2)/(n-1));
    standardized=centered./scale;
    rho=(accumarray(center,zscore(target))'*standardized)/(n-1);
    kept(:,pc)=abs(rho(:))<.01;
    Z=standardized(:,kept(:,pc));assert(size(Z,2)>=k,'Too few retained center draws for endpoint');
    % Compute candidate OLS slopes in blocks without storing N x all draws.
    beta=zeros(size(Z,2),size(Y,2));
    for first=1:128:size(Z,2)
        ix=first:min(first+127,size(Z,2));s=Z(center,ix);s=s-Q*(Q'*s);
        den=sum(s.^2,1)';assert(all(den>1e-10),'Center score is collinear with nuisance design');
        beta(ix,:)=(s'*Y)./den;
    end
    for r=1:size(Y,2)
        features=Z.*beta(:,r)';features=features-counts'*features/n;
        [U,S,~]=svd(sqrt(counts).*features,'econ');v=diag(S);
        assert(numel(v)>=k && v(k)>max(size(features))*eps(v(1)),'Endpoint exceeds artifact rank');
        basis=U(center,1:k)./sqrt(counts(center));
        if nargin>=6
            [P,~]=qr([ones(n,1),protection(:,r,pc)],0);
            basis=basis-P*(P'*basis);
            [A,S,~]=svd(basis,'econ');sv=diag(S);
            A=A(:,sv>max(size(basis))*eps(max(sv)));
            clean(:,r,pc)=zscore(Y(:,r)-A*(A'*Y(:,r)));
        else
            clean(:,r,pc)=zscore(Y(:,r)-basis*(basis'*Y(:,r)));
        end
    end
end
end
