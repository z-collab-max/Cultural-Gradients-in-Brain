function run_china_example(output_dir,n_boot,n_rois)
root=fileparts(fileparts(mfilename('fullpath')));addpath(fileparts(mfilename('fullpath')));
if nargin<1,output_dir=fullfile(root,'outputs','china');end
if nargin<2,n_boot=0;end
if nargin<3,n_rois=3;end
if ~isfolder(output_dir),mkdir(output_dir);end
[T,cortex,behavior,task_names]=make_synthetic_china(20261103,n_rois);
save(fullfile(output_dir,'synthetic_inputs.mat'),'T','cortex','behavior','task_names');
T=prepare_china_table(T);[scores,pca_metadata]=behavior_pcs(behavior);
mixed=fit_china_mixed(T,cortex);relocation=fit_china_relocation(T,cortex);
mediation=fit_china_mediation(T,cortex,scores,n_boot);
behavior_order=best_task_order(mediation.behavior_effect,task_names);
mediation_order=cell(1,2);
for pc=1:2
 if n_boot>1,effects=mediation.task_indirect_z(:,:,pc);else,effects=mediation.task_indirect(:,:,pc);end
 mediation_order{pc}=best_task_order(effects,task_names);
end
roi_names=cellstr(compose('ROI%d',(1:size(cortex,2))'));
save(fullfile(output_dir,'analysis_results.mat'),'mixed','relocation','mediation','behavior_order','mediation_order','pca_metadata','roi_names');
export_china_source(output_dir);plot_china_example(output_dir);
fprintf('CHINA_EXAMPLE_PASS %s\n',output_dir);
end
