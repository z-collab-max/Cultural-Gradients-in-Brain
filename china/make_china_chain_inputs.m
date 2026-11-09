function make_china_chain_inputs(folder)
% Write synthetic keyed files for the shared preprocessing example.
assert(~isfolder(folder),'Choose a new input directory');mkdir(folder);
[T,cortex]=make_synthetic_china(37,6);n=height(T);
id=compose('S%04d',(1:n)');scanner=randi(8,n,1);center_province=scanner;
valid=true(n,1);valid(1:75:end)=false;
participants=table(id,T.province_birth,T.city_birth,T.province_destination,T.relocator,scanner,center_province,valid,...
    'VariableNames',{'participant_id','province_birth','city_birth','province_destination','relocator','scanner','center_province','valid'});
writetable(participants,fullfile(folder,'participants.csv'));
covariates=[table(id,'VariableNames',{'participant_id'}),T(:,china_covariates)];
writetable(covariates(randperm(n),:),fullfile(folder,'covariates.csv'));
[province,first]=unique(T.province_birth);culture=table(province,T.culture_pc1(first),T.culture_pc2(first),...
    'VariableNames',{'province','pc1','pc2'});writetable(culture,fullfile(folder,'culture.csv'));
effects=randn(8,size(cortex,2))*.25;cortex=cortex+effects(scanner,:);
order=randperm(n);cortex=cortex(order,:);participant_id=cellstr(id(order));roi_names=cellstr(compose('ROI%d',1:size(cortex,2)));
save(fullfile(folder,'cortex.mat'),'participant_id','roi_names','cortex','-v7');
end
