function run_china_chain(input_dir,output_dir,outer_draws,inner_draws,seed,protect_destination)
% Keyed individual data, parametric ComBat, center-artifact endpoint and model refits.
if nargin<3,outer_draws=10000;end
if nargin<4,inner_draws=10000;end
if nargin<5,seed=42;end
if nargin<6,protect_destination=true;end
validateattributes(protect_destination,{'logical'},{'scalar'});
validateattributes(outer_draws,{'numeric'},{'integer','positive','scalar'});
validateattributes(inner_draws,{'numeric'},{'integer','positive','scalar'});
assert(~isfolder(output_dir)&&~isfile(output_dir),'Choose a new run directory');
participants=read_keys(fullfile(input_dir,'participants.csv'),...
    {'participant_id','province_birth','city_birth','province_destination','scanner','center_province'});
covariates=read_keys(fullfile(input_dir,'covariates.csv'),{'participant_id'});
culture=read_keys(fullfile(input_dir,'culture.csv'),{'province'});
raw=read_cortex(input_dir);
ids=string(raw.participant_id(:));roi_names=cellstr(string(raw.roi_names(:)));
assert(size(raw.cortex,1)==numel(ids)&&size(raw.cortex,2)==numel(roi_names));
assert(numel(unique(roi_names))==numel(roi_names)&&numel(roi_names)>=2);
assert(all(ismember(china_covariates,covariates.Properties.VariableNames)));
assert(all(ismember(participants.valid,[0 1])));T=participants(logical(participants.valid),:);
ci=row_map(T.participant_id,covariates.participant_id);yi=row_map(T.participant_id,ids);
T=[T,covariates(ci,china_covariates)];Y=double(raw.cortex(yi,:));
assert(all(isfinite(Y),'all')&&all(std(Y)>0),'Incomplete or constant cortical data');
birth=row_map(T.province_birth,culture.province,false);destination=row_map(T.province_destination,culture.province,false);
[centers,~,center]=unique(T.center_province);center_label=row_map(centers,culture.province);
[~,~,batch]=unique(T.scanner);assert(numel(unique(batch))>1&&all(accumarray(batch,1)>=2));
assert(~any(ismissing(T(:,{'scanner','city_birth','province_birth','province_destination','center_province'})),'all'));
scores=culture{:,{'pc1','pc2'}};assert(all(isfinite(scores),'all'));
% Outer rows shuffle paired culture scores across the analysis provinces only.
analysis_units=unique(birth);assert(all(ismember([destination;center_label],analysis_units)));
stream=RandStream('mt19937ar','Seed',seed);inner_orders=zeros(numel(centers),inner_draws);
for j=1:inner_draws,inner_orders(:,j)=randperm(stream,numel(centers))';end
orders=zeros(outer_draws,numel(analysis_units));
for j=1:outer_draws,orders(j,:)=randperm(stream,numel(analysis_units));end
R=size(Y,2);null_beta=zeros(outer_draws,2,R);null_relocation=zeros(outer_draws,R,4,2);
null_groups=zeros(outer_draws,R,2,2);accepted_counts=zeros(outer_draws+1,2);
for b=0:outer_draws
    candidate=scores;
    if b>0,candidate(analysis_units,:)=scores(analysis_units(orders(b,:)),:);end
    for pc=1:2
        T.(sprintf('culture_pc%d',pc))=candidate(birth,pc);
        T.(sprintf('destination_pc%d',pc))=candidate(destination,pc);
    end
    D=prepare_china_table(T);
    protected=[D{:,china_covariates},D.culture_pc1,D.culture_pc2];
    if protect_destination,protected=[protected,D.destination_pc1,D.destination_pc2];end
    harmonized=combat(Y',batch,protected,1)';assert(all(isfinite(harmonized),'all'));
    [clean,kept]=center_artifact_endpoint(harmonized,D,center,candidate(center_label,:),inner_orders);
    mixed=fit_china_mixed(D,clean);reloc=fit_china_relocation(D,clean);
    values=[mixed.beta(:);mixed.se(:);mixed.p(:);reloc.beta(:);reloc.group_beta(:)];
    assert(all(isfinite(values)),'Nonfinite fitted statistic');
    accepted_counts(b+1,:)=sum(kept,1);
    if b==0
        observed_mixed=mixed;observed_relocation=reloc;observed_kept=kept;
        observed_clean=clean;observed_harmonized=harmonized;
        cluster_influence=zeros(2,numel(analysis_units),R);
        for pc=1:2
            webb=prepare_webb_influence(T,clean(:,:,pc),pc);
            assert(max(abs(webb.point-mixed.beta(:,pc)))<1e-5);
            cluster_influence(pc,:,:)=webb.influence;
        end
    else
        null_beta(b,:,:)=mixed.beta';null_relocation(b,:,:,:)=reloc.beta;
        null_groups(b,:,:,:)=reloc.group_beta;
    end
end
parent=fileparts(output_dir);if isempty(parent),parent=pwd;end
if ~isfolder(parent),mkdir(parent);end
stage=tempname(parent);mkdir(stage);
raw_mm=Y;covariates=T{:,china_covariates};province=cellstr(T.province_birth);
scanner=cellstr(T.scanner);birth_scores=scores(birth,:);participant_id=cellstr(T.participant_id);
city=cellstr(string(T.city_birth));
save(fullfile(stage,'fig2_individual_inputs.mat'),'raw_mm','covariates','province','city','scanner','birth_scores','participant_id','roi_names','-v7');
beta=observed_mixed.beta';se=observed_mixed.se';p_model=observed_mixed.p';
p_permutation=reshape(mean(abs(null_beta)>=reshape(abs(beta),1,2,R),1),2,R);
model_labels={'PC1','PC2'};
save(fullfile(stage,'main_maps.mat'),'beta','se','p_model','p_permutation','null_beta','cluster_influence','model_labels','roi_names','-v7');
mover_maps=reshape(observed_relocation.group_beta(:,2,:),R,2)';
stayer_maps=reshape(observed_relocation.group_beta(:,1,:),R,2)';
deculturation_maps=reshape(observed_relocation.beta(:,3,:),R,2)';
acculturation_maps=reshape(observed_relocation.beta(:,4,:),R,2)';
relocation_maps=reshape(observed_relocation.beta(:,2,:),R,2)';
save(fullfile(stage,'commonality.mat'),'mover_maps','stayer_maps','deculturation_maps','acculturation_maps','relocation_maps','-v7');
save(fullfile(stage,'relocation_results.mat'),'observed_relocation','null_relocation','null_groups','roi_names','-v7.3');
participant_id=T.participant_id;province_order=culture.province(analysis_units);endpoint=numel(centers)-1;
save(fullfile(stage,'preprocessing.mat'),'participant_id','roi_names','centers','province_order','endpoint',...
    'orders','inner_orders','observed_kept','accepted_counts','observed_harmonized','observed_clean','seed','protect_destination','-v7.3');
assert(~isfolder(output_dir));[ok,msg]=movefile(stage,output_dir);assert(ok,msg);
end

function t=read_keys(file,names)
opts=detectImportOptions(file,'TextType','string');opts=setvartype(opts,names,'string');t=readtable(file,opts);
end

function raw=read_cortex(folder)
file=fullfile(folder,'cortex.mat');vars=who('-file',file);
if all(ismember({'participant_id','cortex','roi_names'},vars))
    raw=load(file,'participant_id','cortex','roi_names');return
end
assert(ismember('brain_data',vars),'Missing cortical matrix or brain_data structure');
s=load(file,'brain_data');brain=s.brain_data;
names=readlines(fullfile(folder,'roi_names.txt'));names=strip(names);names(names=="")=[];
assert(~isempty(names)&&all(isfield(brain,cellstr(names))),'Missing cortical fields');
ids=brain.MRIID(:);Y=zeros(numel(ids),numel(names));
for r=1:numel(names),Y(:,r)=double(brain.(names(r))(:));end
raw=struct('participant_id',ids,'cortex',Y,'roi_names',names);
end

function ix=row_map(keys,ids,unique_keys)
if nargin<3,unique_keys=true;end
assert(~any(ismissing(ids))&&numel(unique(ids))==numel(ids),'Duplicate or missing input ID');
if unique_keys,assert(numel(unique(keys))==numel(keys),'Duplicate participant ID');end
[ok,ix]=ismember(keys,ids);assert(all(ok),'An eligible ID is missing from an input table');
end

