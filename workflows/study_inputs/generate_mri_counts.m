function generate_mri_counts(brain_file,out_dir)
if ~isfolder(out_dir),mkdir(out_dir);end

%%


dem = 2;

%% get preparation
string_dem = {'area','thickness','volume'};

load(brain_file)% calculate total area
all_vars = fieldnames(brain_data);
area_var = all_vars(contains(all_vars,'area') & ~contains(all_vars,'500_aparc_'));
areadata = nan(length(brain_data.age),length(area_var));
for arvv = 1:length(area_var)
    areadata(:,arvv) = brain_data.(area_var{arvv});
end
brain_data.Total_area = sum(areadata,2);
controlvars = {'sex','age','edu','Gpc1','Gpc2','Gpc3','Gpc4','Gpc5','Gpc6','Gpc7','Gpc8','Gpc9','Gpc10','Incomeindex21','Educationalindex21',...
   'Incomeindex21_recruit','Educationalindex21_recruit'}; %,'Incomeindex21_recruit','Educationalindex21_recruit'

% do not change the 'sex' order!!
naturalvars = {'C_sunshine_18y','C_humidity_18y','C_temperature_18y',...
    'population_density_2015_b','Green_Covered_Area_prop_b'};
Culnames_birth = {'pc7_1','pc7_2','pc7_3'};
Culnames_recruit = {'pc7_1_recruit','pc7_2_recruit'};

controlvars = [controlvars naturalvars];
all_vars = fieldnames(brain_data);
if dem == 2
    brainidx = (contains(all_vars,string_dem{dem}) & ~contains(all_vars,{'500_aparc_','unknown','WhiteSurfArea'}));
elseif dem == 3
    brainidx = (contains(all_vars,string_dem{dem}) & ~contains(all_vars,{'500_aparc_','unknown','WhiteSurfArea'}))...
        & ~((contains(all_vars,'TIV') & ~(contains(all_vars,'eTIV'))));
    controlvars = [controlvars 'TIV'];
elseif dem == 1
    brainidx = (contains(all_vars,string_dem{dem}) & ~contains(all_vars,{'500_aparc_','unknown','WhiteSurfArea'}))...
        & ~(contains(all_vars,'Total'));
    controlvars = [controlvars 'Total_area'];
end
brain_var = all_vars(brainidx);
controlM = nan(length(brain_data.age),length(controlvars));
for cvv = 1:length(controlvars)
    controlM(:,cvv) = brain_data.(controlvars{cvv})(:);
end
brainM = nan(length(brain_data.age),length(brain_var));
for bvv = 1:length(brain_var)
    brainM(:,bvv) = brain_data.(brain_var{bvv});
end
predM_birth = nan(length(brain_data.age),length(Culnames_birth));
for hvv = 1:length(Culnames_birth)
    predM_birth(:,hvv) = brain_data.(Culnames_birth{hvv});
end

predM_recruit = nan(length(brain_data.age),length(Culnames_recruit));
for hvv = 1:length(Culnames_recruit)
    predM_recruit(:,hvv) = brain_data.(Culnames_recruit{hvv});
end

valididx = ~isnan(sum(controlM,2)) & ~isnan(sum(predM_birth,2)) & ~isnan(sum(predM_recruit,2)) & ~isnan(sum(brainM,2)) & ~isnan(brain_data.centerID_rerange);%...
center_code = brain_data.centerID_rerange(valididx);


province_names = cellstr(brain_data.province_birth(valididx));
[names,~,groups]=unique(province_names);
counts=accumarray(groups,1);
T=table(string(names),counts,'VariableNames',{'Province','Subject_Count'});
T=sortrows(T,'Subject_Count','descend');
writetable(T,fullfile(out_dir,'Figure1b_province_counts.csv'),'Encoding','UTF-8');
city_names=cellstr(brain_data.city_birth(valididx));
[names,~,groups]=unique(city_names);counts=accumarray(groups,1);
C=table(string(names),counts,log(counts),'VariableNames',{'City','Subject_Count','Log_Subject_Count'});
writetable(sortrows(C,'Subject_Count','descend'),fullfile(out_dir,'MRI_cities.csv'),'Encoding','UTF-8');
summary=struct('n',sum(valididx),'province_count',height(T),'excluded_variables',{{'Droughts_b','Floods_b','Terrain_Ruggedness_b','Agricultural_Suitability_b'}});
fid=fopen(fullfile(out_dir,'MRI_valid_summary.json'),'w');fprintf(fid,'%s',jsonencode(summary));fclose(fid);
fprintf('MRI valid n=%d, provinces=%d\n',sum(valididx),height(T));

end
