function names=china_covariates(T)
names=[{'age','sex','education'},arrayfun(@(k)sprintf('genetic_pc%d',k),1:10,'UniformOutput',false),...
 {'father_education','mother_education','intracranial_volume','total_surface_area'}];
if nargin==0
 macro=[arrayfun(@(k)sprintf('city_macro_pc%d',k),1:9,'UniformOutput',false),...
 arrayfun(@(k)sprintf('province_macro_pc%d',k),1:2,'UniformOutput',false)];
else
 if istable(T),fields=T.Properties.VariableNames;else,fields=fieldnames(T)';end
 macro={};
 for prefix={'city_macro_pc','province_macro_pc'}
  matched=fields(~cellfun('isempty',regexp(fields,['^' prefix{1} '[1-9][0-9]*$'],'once')));
  assert(~isempty(matched),'Missing macro PCs; run prepare_macro_pcs first');
  index=cellfun(@(v)str2double(extractAfter(v,prefix{1})),matched);[index,order]=sort(index);
  assert(isequal(index,1:numel(index)),'PC indices must be consecutive');macro=[macro matched(order)];
 end
 assert(all(ismember(names,fields)),'Missing individual covariates');
 assert(~any(ismember({'absolute_latitude','hdi_income'},fields)),'Run prepare_macro_pcs before model fitting');
end
names=[names macro];
end
