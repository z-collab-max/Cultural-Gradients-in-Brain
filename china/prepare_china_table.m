function T=prepare_china_table(T)
required=[china_covariates,{'province_birth','city_birth','culture_pc1','culture_pc2',...
 'relocator','destination_pc1','destination_pc2'}];
assert(all(ismember(required,T.Properties.VariableNames)),'Required input fields are missing.');
assert(~any(ismissing(T(:,required)),'all'),'Input table contains missing values.');
numeric=[china_covariates,{'culture_pc1','culture_pc2','relocator','destination_pc1','destination_pc2'}];
assert(all(isfinite(T{:,numeric}),'all'),'Numeric input fields must be finite.');
assert(all(ismember(T.relocator,[0 1])),'relocator must be 0 or 1.');
T.province_city=categorical(string(T.province_birth)+"__"+string(T.city_birth));
T.province_birth=categorical(T.province_birth);
names=china_covariates;
for j=1:numel(names)
 x=double(T.(names{j}));s=std(x,'omitnan');assert(s>0);
 T.(names{j})=(x-mean(x,'omitnan'))/s;
end
for pc=1:2
 name=sprintf('culture_pc%d',pc);destination=sprintf('destination_pc%d',pc);
 mu=mean(T.(name),'omitnan');s=std(T.(name),'omitnan');assert(s>0);
 T.(destination)=(T.(destination)-mu)/s;T.(name)=(T.(name)-mu)/s;
 delta=T.(destination)-T.(name);delta(T.relocator==0)=0;
 T.(sprintf('culture_delta%d',pc))=delta;
end
end
