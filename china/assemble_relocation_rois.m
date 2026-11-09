function assemble_relocation_rois(S,out,B)
% Aligned, display-oriented birth/destination scores, harmonized thickness and gradients.
if nargin<3,B=10000;end
validateattributes(B,{'numeric'},{'scalar','integer','>=',2});
names=cellstr(string(S.roi_names(:)));N=numel(S.participant_id);R=numel(names);
assert(numel(unique(string(S.participant_id)))==N && numel(unique(names))==R);
assert(isequal(size(S.thickness),[N R]) && size(S.covariates,1)==N && isequal(size(S.birth),[N 2]) && isequal(size(S.destination),[N 2]));
assert(isequal(size(S.beta),[R 2]) && isequal(size(S.gradients),[R 2]));
X=[ones(N,1),S.covariates];assert(all(isfinite([X,S.thickness,S.birth,S.destination]),'all') && all(isfinite([S.gradients,S.beta]),'all'));
assert(all(ismember(S.mover,[0 1]))&&numel(S.mover)==N);
Y=zscore(S.thickness-X*(X\S.thickness));assert(all(isfinite(Y),'all'));
groups={'HH','LL','HL','LH','MHH','MLL'};mover=logical(S.mover(:));group=zeros(N,2);network_means=nan(N,2,2);
posidx=false(R,2);negidx=posidx;thresholds=zeros(2,4);roi_effect=zeros(R,2);
for pc=1:2
    grad=S.gradients(:,pc);beta=S.beta(:,pc);thresholds(pc,:)=[prctile(grad,[33 66]),prctile(beta,[10 90])];
    t=thresholds(pc,:);posidx(:,pc)=grad>t(2)&beta>t(4);negidx(:,pc)=grad<t(1)&beta<t(3);
    assert(any(posidx(:,pc))&&any(negidx(:,pc)),'Empty selection for PC%d; check map order/orientation',pc);
    birth=S.birth(:,pc);dest=S.destination(:,pc);cut=median(birth);
    if pc==1,high=birth>=cut;arrive=dest>=cut;else,high=birth>cut;arrive=dest>cut;end
    mask=[~mover&high&arrive,~mover&~high&~arrive,mover&high&~arrive,mover&~high&arrive,mover&high&arrive,mover&~high&~arrive];
    assert(all(sum(mask,2)<=1)&&all(sum(mask)>2),'Each of the six groups needs three observations');
    group(:,pc)=mask*(1:6)';means=[mean(Y(:,posidx(:,pc)),2),mean(Y(:,negidx(:,pc)),2)];network_means(:,:,pc)=means;
    g=cell(1,6);for k=1:6,g{k}=means(mask(:,k),:);end
    ns=sum(mask);n=min(ns(3:4));rng(20261027+10*pc+1,'twister');draw=zeros(B,2);
    for j=1:B
        ih=randsample(ns(1),n,true);il=randsample(ns(2),n,true);
        draw(j,:)=pooled_d(g{1}(ih,:),g{2}(il,:));
    end
    cross=pooled_d(g{3},g{4});same=pooled_d(g{5},g{6});
    d=struct;d.treatcohendMP=cross(1);d.treatcohendMN=cross(2);d.treatcohend_MHH_MLL_P=same(1);d.treatcohend_MHH_MLL_N=same(2);
    d.cohendpool_control_posM=draw(:,1);d.cohendpool_control_negM=draw(:,2);d.diff_dist=draw(:,1)-draw(:,2);
    d.diff_treat_cross=cross(1)-cross(2);d.diff_treat_same=same(1)-same(2);obs=[d.diff_treat_cross,d.diff_treat_same];
    p=min(1,2*min(mean(d.diff_dist<=obs),mean(d.diff_dist>=obs)));d.p_val_cross=p(1);d.p_val_same=p(2);
    d.group_n=ns;d.birth_median=cut;d.sample_size=n;
    rng(20261104+pc,'twister');
    for side={'pos','neg'}
        tag=side{1};if strcmp(tag,'pos'),ids=find(posidx(:,pc));direction='descend';else,ids=find(negidx(:,pc));direction='ascend';end
        [stable,ci,transition,dp]=roi_intervals(Y(:,ids),mask,B);
        [~,order]=sort(stable,direction);[~,mean_order]=sort(dp(:,5)-dp(:,1),direction);
        d.(['idx_' tag '_sorted'])=order;d.(['mean_stable_' tag])=stable;d.(['ci_stable_' tag])=ci;d.(['obs_d_' tag])=transition;
        d.(['dp_data_' tag])=dp(mean_order,:);d.(['dp_names_' tag])=names(ids(mean_order));d.(['effect_names_' tag])=names(ids(order));
    end
    Fig4_ED9_Data(pc)=d;roi_effect(:,pc)=pooled_d(Y(mask(:,1),:),Y(mask(:,2),:))';
end
% ED11a uses a single random stream across the two PCs.
rng(20261026,'twister');panels=cell(1,2);
for pc=1:2
    values=network_means(:,1,pc)-network_means(:,2,pc);g=cell(1,6);
    for j=1:6,g{j}=values(group(:,pc)==j);end
    n=min(numel(g{3}),numel(g{4}));H=zeros(B,1);L=H;
    for j=1:B,H(j)=mean(randsample(g{1},n,true));L(j)=mean(randsample(g{2},n,true));end
    E.null_H=H;E.null_L=L;E.null_difference=H-L;E.groups=groups;E.group_n=cellfun(@numel,g);E.means=cellfun(@mean,g);E.sample_size=n;
    E.cross=E.means(3)-E.means(4);E.same=E.means(5)-E.means(6);tail=@(v)max(1/B,min(mean(H-L>v),mean(H-L<v)));
    E.p_cross=tail(E.cross);E.p_same=tail(E.same);ED(pc)=E;panels{pc}=E;panels{pc}.folder='ExtendedDataFigure11';panels{pc}.name=sprintf('Network_PC%d',pc);panels{pc}.kind='paired_histogram';
end
if ~isfolder(out),mkdir(out);end
colorpair={'#f4a582','#E6A094','#d6604d';'#92c5de','#8DA1C6','#4163a0'};
birthplace_beta=S.beta;canonical_gradient=S.gradients;
save(fullfile(out,'roi_source.mat'),'Fig4_ED9_Data','colorpair','roi_effect','names','posidx','negidx','thresholds','birthplace_beta','canonical_gradient','B','-v7');
save(fullfile(out,'network_source.mat'),'ED','panels','B','-v7');
participant_id=S.participant_id;save(fullfile(out,'individual_roi_means.mat'),'participant_id','group','network_means','groups','-v7');
end

function d=pooled_d(a,b)
na=size(a,1);nb=size(b,1);d=(mean(a)-mean(b))./sqrt(((na-1)*var(a)+(nb-1)*var(b))/(na+nb-2));
end

function [stable,ci,transition,dp]=roi_intervals(Y,mask,B)
nr=size(Y,2);stable=zeros(nr,1);ci=zeros(nr,2);transition=stable;dp=zeros(nr,8);
n=round(mean(sum(mask(:,3:4))));
for r=1:nr
    hh=Y(mask(:,1),r);ll=Y(mask(:,2),r);hl=Y(mask(:,3),r);lh=Y(mask(:,4),r);
    transition(r)=(mean(hl)-mean(lh))/sqrt((var(hl)+var(lh))/2);
    ds=zeros(B,1);mh=ds;ml=ds;
    for start=1:128:B
        ix=start:min(start+127,B);q=numel(ix);high=hh(randi(numel(hh),n,q));low=ll(randi(numel(ll),n,q));
        ds(ix)=(mean(high,1)-mean(low,1))./sqrt((var(high,0,1)+var(low,0,1))/2);
        mh(ix)=mean(hh(randi(numel(hh),numel(hh),q)),1);ml(ix)=mean(ll(randi(numel(ll),numel(ll),q)),1);
    end
    stable(r)=mean(ds);ci(r,:)=prctile(ds,[2.5 97.5]);
    dp(r,:)=[mean(ll),mean(lh),prctile(ml,[2.5 97.5]),mean(hh),mean(hl),prctile(mh,[2.5 97.5])];
end
end
