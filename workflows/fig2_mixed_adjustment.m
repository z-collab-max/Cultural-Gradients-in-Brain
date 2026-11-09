function fig2_mixed_adjustment(input_file,output_file)
% Adjust parcel values using the shared China REML fixed covariate estimates.
here=fileparts(mfilename('fullpath'));addpath(fullfile(here,'study_inputs'),fullfile(here,'..','china'));
S=load(input_file);names=china_covariates;b=struct;
assert(size(S.covariates,2)==28 && size(S.raw,2)==numel(S.pc_index));
assert(all(isfinite([S.raw S.covariates S.scores S.groups S.scanner]),'all'));
for k=1:28,b.(names{k})=S.covariates(:,k);end
b.province_birth=S.groups(:,1);b.province_city=S.groups(:,2);T=study_china_table(b);
assert(numel(categories(T.province_birth))>=2 && numel(categories(T.province_city))>=2);
[~,~,scanner]=unique(S.scanner);adjusted_mm=zeros(size(S.raw));beta=zeros(1,size(S.raw,2));p=beta;
formula=china_formula('response',{'culture'});
for r=1:size(S.raw,2)
 raw=S.raw(:,r);center_mean=accumarray(scanner,raw,[],@mean);y=raw-center_mean(scanner);
 assert(std(y)>0 && std(raw)>0);
 D=T;D.response=zscore(y);D.culture=zscore(S.scores(:,S.pc_index(r)));
 assert(rank([ones(height(D),1),D{:,names},D.culture])==30,'Rank-deficient fixed design');
 m=fitlme(D,formula,'FitMethod','REML');
 [found,indices]=ismember(names,m.Coefficients.Name);assert(all(found));
 adjusted_z=D.response-D{:,names}*m.Coefficients.Estimate(indices);
 adjusted_mm(:,r)=mean(raw)+std(raw)*adjusted_z;
 j=strcmp(m.Coefficients.Name,'culture');beta(r)=m.Coefficients.Estimate(j);p(r)=m.Coefficients.pValue(j);
end
save(output_file,'adjusted_mm','beta','p','formula','-v7');
end
