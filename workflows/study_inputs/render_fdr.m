function render_fdr(main_maps,out)
% Benjamini-Hochberg correction within each cultural PC across tested ROIs.
here=fileparts(mfilename('fullpath'));addpath(fullfile(here,'..','..','figures'));
if ismember('AllRes',who('-file',main_maps))
    raw=load(main_maps,'AllRes');t=raw.AllRes.Expcs.TValue(:,1:2)';roi_names=raw.AllRes.Info.BrainVars;
    df=29;
else
    s=load(main_maps,'beta','se','roi_names','cluster_influence');t=s.beta./s.se;roi_names=s.roi_names;
    df=size(s.cluster_influence,2)-1;
end
if ~isfolder(out),mkdir(out);end
for j=1:size(t,1)
    values=t(j,:)'*(3-2*j);p=2*tcdf(-abs(values),df);valid=isfinite(p);
    assert(all(p(valid)>=0 & p(valid)<=1));idx=find(valid);[sorted,i]=sort(p(valid));
    m=numel(idx);q=nan(size(p));q(idx(i))=min(1,flip(cummin(flip(sorted.*m./(1:m)'))));
    significant=q<0.05;color_limits=[-1 1]*(6-j);
    stem=sprintf('PC%d_FDR',j);file=fullfile(out,[stem '.mat']);
    save(file,'roi_names','values','p','q','significant','color_limits');render_brain_map(file,fullfile(out,stem));
end
end


