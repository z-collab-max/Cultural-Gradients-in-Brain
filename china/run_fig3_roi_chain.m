function run_fig3_roi_chain(input_dir,results,spatial_file,out,B)
% Build Fig.3 and ED11 sources from one shared China run and two aligned gradients.
if nargin<5,B=10000;end
P=load(fullfile(results,'preprocessing.mat'),'observed_harmonized','participant_id','roi_names');
R=load(fullfile(results,'relocation_results.mat'),'observed_relocation','roi_names');
G=load(spatial_file,'roi_names','gradients','pc_sign');names=string(P.roi_names(:));
assert(isequal(names,string(R.roi_names(:)))&&isequal(names,string(G.roi_names(:))),'ROI order mismatch');
signs=double(G.pc_sign(:))';assert(numel(signs)==2&&all(ismember(signs,[-1 1])));
T=read_ids(fullfile(input_dir,'participants.csv'),{'participant_id','province_birth','province_destination'});
C=read_ids(fullfile(input_dir,'covariates.csv'),{'participant_id'});U=read_ids(fullfile(input_dir,'culture.csv'),{'province'});
ids=string(P.participant_id(:));T=T(locate(ids,T.participant_id),:);C=C(locate(ids,C.participant_id),:);
birth=locate(T.province_birth,U.province);dest=locate(T.province_destination,U.province);
S.participant_id=ids;S.thickness=P.observed_harmonized;S.roi_names=names;S.covariates=C{:,china_covariates};S.mover=T.relocator;
S.birth=U{birth,{'pc1','pc2'}}.*signs;S.destination=U{dest,{'pc1','pc2'}}.*signs;
S.beta=reshape(R.observed_relocation.beta(:,1,:),[],2).*signs;S.gradients=G.gradients;
assemble_relocation_rois(S,out,B);
end

function t=read_ids(file,names)
o=detectImportOptions(file,'TextType','string');o=setvartype(o,names,'string');t=readtable(file,o);
end

function i=locate(keys,values)
assert(~any(ismissing(values))&&numel(unique(values))==numel(values),'Duplicate or missing input key');
[ok,i]=ismember(keys,values);assert(all(ok),'Missing input key');
end
