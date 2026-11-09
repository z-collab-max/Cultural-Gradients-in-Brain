function render_study_outputs(input_root,output_root)
% Render each panel recipe and named surface input beneath a run directory.
here=fileparts(mfilename('fullpath'));addpath(here);
files=dir(fullfile(input_root,'**','panels.mat'));
for j=1:numel(files)
 relative=erase(files(j).folder,[char(java.io.File(input_root).getCanonicalPath) filesep]);
 if strcmp(relative,char(java.io.File(input_root).getCanonicalPath)),relative='';end
 render_panels(fullfile(files(j).folder,files(j).name),fullfile(output_root,relative));
end
files=dir(fullfile(input_root,'**','*surface.mat'));
for j=1:numel(files)
 relative=erase(files(j).folder,[char(java.io.File(input_root).getCanonicalPath) filesep]);
 if strcmp(relative,char(java.io.File(input_root).getCanonicalPath)),relative='';end
 [~,name]=fileparts(files(j).name);render_brain_map(fullfile(files(j).folder,files(j).name),fullfile(output_root,relative,name));
end
end
