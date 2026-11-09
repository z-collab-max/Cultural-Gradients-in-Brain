function study_china_batch(input_file,output_file)
S=load(input_file);addpath(fullfile(fileparts(mfilename('fullpath')),'..','..','china'));
names=china_covariates;b=struct;
for k=1:28,b.(names{k})=S.X(:,k);end
b.province_birth=S.groups(:,1);b.province_city=S.groups(:,2);
T=study_china_table(b);result=study_china_fit(T,S.Y,S.predictors);
save(output_file,'-struct','result','-v7');
end
