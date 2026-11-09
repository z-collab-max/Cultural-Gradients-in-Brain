function generate_brain_icc(brain_file,output_dir)
if ~isfolder(output_dir),mkdir(output_dir);end
load(brain_file,'brain_data');
addpath(fileparts(mfilename('fullpath')));
addpath(fullfile(fileparts(mfilename('fullpath')),'..','..','china'));
ExtFig1_ICC_Data = struct();
string_dem = {'area','thickness','volume'};

for dem = 1:3
    fprintf('\n=========================================\n');
    fprintf('Processing morphological indicators: %s ...\n', string_dem{dem});
    
    all_vars=fieldnames(brain_data);
    controlvars=china_covariates;
    Culnames_birth={'pc7_1','pc7_2'};Culnames_recruit={};
    % 2. Select Brain Features based on modality
    if dem == 2
        brainidx = (contains(all_vars,string_dem{dem}) & ~contains(all_vars,{'500_aparc_','unknown','WhiteSurfArea'}));
    elseif dem == 3
        brainidx = (contains(all_vars,string_dem{dem}) & ~contains(all_vars,{'500_aparc_','unknown','WhiteSurfArea'}))...
            & ~((contains(all_vars,'TIV') & ~(contains(all_vars,'eTIV'))));

    elseif dem == 1
        brainidx = (contains(all_vars,string_dem{dem}) & ~contains(all_vars,{'500_aparc_','unknown','WhiteSurfArea'}))...
            & ~(contains(all_vars,{'Total','total_surface_area'}));

    end
    brain_var = all_vars(brainidx & ~ismember(all_vars,china_covariates));
    
    % 3. Build Matrices
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
    
    valididx = ~isnan(sum(controlM,2)) & ~isnan(sum(predM_birth,2)) & ~isnan(sum(predM_recruit,2)) & ~isnan(sum(brainM,2)) & ~isnan(brain_data.centerID_rerange);
    center_code = brain_data.centerID_rerange(valididx);
    n_cul = length(Culnames_birth);
    
    [~,~,site]=unique(center_code);data=zeros(sum(valididx),numel(brain_var));
    for r=1:numel(brain_var)
        y=brainM(valididx,r);means=accumarray(site,y,[],@mean);data(:,r)=zscore(y-means(site));
    end
    T=study_china_table(brain_data,find(valididx));
    data2=zscore(study_china_adjust(T,data,predM_birth(valididx,:)));
    % 6. Calculate ICC1 & ICC3
    fprintf('  - Calculate intraclass correlation coefficient (ICC1 & ICC3)...\n');
    grouplabel = brain_data.province_birth(valididx);
    icc1_vec = nan(size(brain_var));
    icc3_vec = nan(size(brain_var));
    for bvv = 1:length(brain_var)
        [icc1_vec(bvv),icc3_vec(bvv)] = ICC1(data2(:,bvv), grouplabel);
    end
    
    % 7. Save to Structure
    ExtFig1_ICC_Data.(string_dem{dem}).icc1 = icc1_vec;
    ExtFig1_ICC_Data.(string_dem{dem}).icc3 = icc3_vec;
    calculation.(string_dem{dem}).valididx=valididx;
    calculation.(string_dem{dem}).brain_var=brain_var;
    calculation.(string_dem{dem}).controlvars=controlvars;
    calculation.(string_dem{dem}).n=sum(valididx);
    calculation.(string_dem{dem}).total_area=brain_data.Total_area;
    save(fullfile(output_dir,'ICC_calculation.mat'),'calculation','-v7.3');
    
    fprintf('  > % s modal calculation is completed and stored in the structure. \n', string_dem{dem});
end

save(fullfile(output_dir,'ExtFig1_ICC_SourceData.mat'),'ExtFig1_ICC_Data','string_dem');
end
