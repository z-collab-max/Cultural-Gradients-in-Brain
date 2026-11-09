function T=study_china_table(b,rows)
% Canonical 28 covariates and birthplace province/city identifiers.
addpath(fullfile(fileparts(mfilename('fullpath')),'..','..','china'));
if nargin<2,rows=(1:numel(b.age))';end
names=china_covariates;T=table;
for k=1:numel(names)
 assert(isfield(b,names{k}),'Missing canonical field: %s',names{k});
 v=b.(names{k});T.(names{k})=reshape(double(v(rows)),[],1);
end
assert(all(isfinite(T{:,:}),'all'),'Covariates must be finite');
assert(all(std(T{:,:})>0),'Covariates must have positive variance');
T{:,:}=zscore(T{:,:});
province=reshape(b.province_birth(rows),[],1);city=reshape(b.province_city(rows),[],1);
assert(~any(ismissing(province)) && ~any(ismissing(city)),'Missing group identifiers');
T.province_birth=categorical(province);
T.province_city=categorical(findgroups(T.province_birth,categorical(city)));
end
