function generate_participant_inputs(brain_file,valid_file,output_dir)
% Export participant identifiers and province labels for residence processing.
if ~isfolder(output_dir),mkdir(output_dir);end
B=load(brain_file,'brain_data');b=B.brain_data;
V=load(valid_file,'valididx_base');
id=double(b.MRIID(:));valid=logical(V.valididx_base(:));lived=double(b.livedplace_32(:));
birth=string(b.province_birth(:));destination=string(b.province_recruit(:));
row=(1:numel(id))';T=table(row,id,valid,lived,birth,destination);
writetable(T,fullfile(output_dir,'participant_inputs.csv'),'Encoding','UTF-8');
end
