function export_china_source(output_dir)
S=load(fullfile(output_dir,'analysis_results.mat'));
source=fullfile(output_dir,'source_data');if ~isfolder(source),mkdir(source);end
roi=(1:size(S.mixed.beta,1))';
for pc=1:2
 beta=S.mixed.beta(:,pc);se=S.mixed.se(:,pc);p=S.mixed.p(:,pc);
 writetable(table(roi,beta,se,p),fullfile(source,sprintf('mixed_PC%d.csv',pc)));
 birthplace=S.relocation.beta(:,1,pc);relocation_shift=S.relocation.beta(:,2,pc);
 deculturation=S.relocation.beta(:,3,pc);acculturation=S.relocation.beta(:,4,pc);
 stayer=S.relocation.group_beta(:,1,pc);relocator=S.relocation.group_beta(:,2,pc);
 difference=stayer-relocator;
 writetable(table(roi,birthplace,relocation_shift,deculturation,acculturation,stayer,relocator,difference),...
  fullfile(source,sprintf('relocation_PC%d.csv',pc)));
 order=S.mediation_order{pc};rsa=order.rsa_map;slope=order.slope;
 writetable(table(roi,rsa,slope),fullfile(source,sprintf('mediation_PC%d.csv',pc)));
 save(fullfile(source,sprintf('mediation_RDM_PC%d.mat',pc)),'order');
end
order=S.behavior_order;rsa=order.rsa_map;slope=order.slope;
writetable(table(roi,rsa,slope),fullfile(source,'brain_behavior.csv'));
save(fullfile(source,'behavior_RDM.mat'),'order');
mediation=S.mediation;behavior_order=S.behavior_order;
save(fullfile(source,'hierarchy_inputs.mat'),'mediation','behavior_order','-v7');
if isfield(S,'roi_names'),roi_names=S.roi_names;save(fullfile(source,'hierarchy_inputs.mat'),'roi_names','-append');end
end
