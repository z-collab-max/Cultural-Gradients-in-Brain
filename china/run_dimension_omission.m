function run_dimension_omission(input_file,output_file)
% Input: T, cortex, cultural_dimensions [province x dimension], province_ids.
S=load(input_file);raw=S.T;Y=S.cortex;C=double(S.cultural_dimensions);
[ok,group]=ismember(string(raw.province_birth),string(S.province_ids));assert(all(ok));
[~,full]=pca(zscore(C));full=full(:,1:2);D=size(C,2);assert(D>=3);
T=raw;for pc=1:2,T.(sprintf('culture_pc%d',pc))=full(group,pc);end
fit=fit_china_mixed(prepare_china_table(T),Y);full_maps=fit.beta';
fold_maps=zeros(D,2,size(Y,2));matching=zeros(D,2);
for omit=1:D
 [~,score]=pca(zscore(C(:,setdiff(1:D,omit))));cross=corr(full,score(:,1:2));
 if abs(cross(1,1))+abs(cross(2,2))>=abs(cross(1,2))+abs(cross(2,1)),order=[1 2];else,order=[2 1];end
 matching(omit,:)=order;
 for pc=1:2
  sign_axis=sign(cross(pc,order(pc)));if sign_axis==0,sign_axis=1;end
  T.(sprintf('culture_pc%d',pc))=score(group,order(pc))*sign_axis;
 end
 fit=fit_china_mixed(prepare_china_table(T),Y);fold_maps(omit,:,:)=fit.beta';
end
if isfield(S,'dimension_labels'),dimension_labels=S.dimension_labels;else,dimension_labels=compose('Dimension%d',1:D);end
parent=fileparts(output_file);if ~isempty(parent)&&~isfolder(parent),mkdir(parent);end
save(output_file,'full_maps','fold_maps','matching','dimension_labels','-v7');
end
