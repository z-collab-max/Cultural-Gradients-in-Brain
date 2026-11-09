function prepare_combat_inputs(brain_file,culture_file,output_file)
% Prepare individual thickness, nuisance variables, and ComBat variants.
here=fileparts(mfilename('fullpath'));

S=load(brain_file,'brain_data');d=S.brain_data;
names=fieldnames(d);roi=names(contains(names,'thickness') & ~contains(names,{'500_aparc_','unknown','WhiteSurfArea'}));
controls=china_covariates;
X=columns(d,controls);Y=columns(d,roi);birth=columns(d,{'pc7_1','pc7_2','pc7_3'});
recruit=columns(d,{'pc7_1_recruit','pc7_2_recruit'});
assert(isnumeric(d.province_birth)&&isnumeric(d.province_city),'Use numeric group identifiers');
valid=all(isfinite([X Y birth recruit double(d.centerID_rerange(:)) double(d.province_birth(:)) double(d.province_city(:))]),2);
X=X(valid,:);Y=Y(valid,:);birth=birth(valid,:);recruit=recruit(valid,:);center=d.centerID_rerange(valid);
assert(all(isfinite([X Y birth recruit]),'all'));
gender=dummyvar(categorical(X(:,2)));B=dummyvar(categorical(center));Y_center=zeros(size(Y));
for j=1:size(Y,2),[~,~,stats]=glmfit(B(:,1:end-1),Y(:,j));Y_center(:,j)=zscore(stats.resid);end
harmonized=struct;
for p=2:3
 model=[gender(:,2) X(:,[1 3:13]) birth(:,1:p) recruit];
 harmonized.(sprintf('PC%d_recruit1',p))=zscore(combat(Y',center,model,1)');
end
harmonized.no_culture=zscore(combat(Y',center,[gender(:,2) X(:,[1 3:13])],1)');
model_design=struct;for k=1:28,model_design.(controls{k})=X(:,k);end
model_design.province_birth=d.province_birth(valid);model_design.province_city=d.province_city(valid);
assert(isnumeric(model_design.province_birth)&&isnumeric(model_design.province_city),'Use numeric group identifiers');
C=load(culture_file,'Culdata');culture=[C.Culdata.pc7_1(:) C.Culdata.pc7_2(:)];
out=fileparts(output_file);if ~isfolder(out),mkdir(out);end
save(output_file,'X','birth','recruit','center','valid','roi','controls','Y_center','harmonized','culture','model_design','-v7.3');
end
function X=columns(d,names)
X=zeros(numel(d.age),numel(names));
for j=1:numel(names),X(:,j)=d.(names{j})(:);end
end
