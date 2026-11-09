function [T,Y,center,center_scores,orders,roi_names]=load_china_run(input_dir,run_dir)
% Align a completed shared run to its keyed participant, covariate and culture tables.
s=load(fullfile(run_dir,'preprocessing.mat'));
T=readtable(fullfile(input_dir,'participants.csv'),'TextType','string');
C=readtable(fullfile(input_dir,'covariates.csv'),'TextType','string');
U=readtable(fullfile(input_dir,'culture.csv'),'TextType','string');
id=string(s.participant_id(:));T=T(align(id,string(T.participant_id)),:);
C=C(align(id,string(C.participant_id)),:);
T=[T,C(:,china_covariates)];
birth=align(string(T.province_birth),string(U.province));
dest=align(string(T.province_destination),string(U.province));
centers=string(s.centers);center=align(string(T.center_province),centers);
unit=align(centers,string(U.province));center_scores=U{unit,{'pc1','pc2'}};
for pc=1:2
    field=sprintf('pc%d',pc);
    T.(sprintf('culture_pc%d',pc))=U.(field)(birth);
    T.(sprintf('destination_pc%d',pc))=U.(field)(dest);
end
T=prepare_china_table(T);Y=s.observed_harmonized;orders=s.inner_orders;
roi_names=s.roi_names;
assert(height(T)==size(Y,1));
end

function indices=align(keys,reference)
assert(numel(unique(reference))==numel(reference)&&~any(ismissing(reference)));
[found,indices]=ismember(keys,reference);assert(all(found),'Unmatched input key');
end
