function prepare_brain_macro_pcs(input_file,output_file)
% Preserve participant order while replacing raw macro fields in brain_data.
assert(~strcmp(input_file,output_file),'Use a separate output file');
S=load(input_file,'brain_data');brain_data=S.brain_data;
raw={'absolute_latitude','terrain_ruggedness','historical_humidity','historical_temperature',...
 'historical_sunshine','historical_precipitation','population_density','green_space',...
 'natural_disaster_frequency_per_capita','hdi_education','hdi_income'};
T=table(brain_data.province_birth(:),brain_data.province_city(:),'VariableNames',{'province_birth','city_birth'});
for j=1:numel(raw),T.(raw{j})=double(brain_data.(raw{j})(:));end
[T,macro_pca]=prepare_macro_pcs(T);brain_data=rmfield(brain_data,raw);
names=[macro_pca.city.names macro_pca.province.names];
for j=1:numel(names),brain_data.(names{j})=T.(names{j});end
parent=fileparts(output_file);if ~isempty(parent)&&~isfolder(parent),mkdir(parent);end
save(output_file,'brain_data','macro_pca','-v7.3');
end
