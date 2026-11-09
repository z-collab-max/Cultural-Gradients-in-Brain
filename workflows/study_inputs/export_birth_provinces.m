function export_birth_provinces(brain_file,out)
s=load(brain_file,'brain_data');b=s.brain_data;
subject_id=string(b.MRIID(:));province=string(b.province_birth(:));
assert(numel(unique(subject_id))==numel(subject_id));
writetable(table(subject_id,province),out,'Encoding','UTF-8');
end
