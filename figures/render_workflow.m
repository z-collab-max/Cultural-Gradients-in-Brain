function render_workflow(source_dir,output_dir)
% Execute a source directory's render_manifest.json.
here=fileparts(mfilename('fullpath'));addpath(here);addpath(fullfile(here,'layouts'));
if nargin<2,output_dir=fullfile(source_dir,'figures');end
jobs=jsondecode(fileread(fullfile(source_dir,'render_manifest.json')));
for k=1:numel(jobs)
 if iscell(jobs),j=jobs{k};else,j=jobs(k);end
 assert(~contains(j.input,'..')&&isempty(regexp(j.input,'[\\/:]','once')));
 input=fullfile(source_dir,j.input);
 switch j.kind
  case 'native',render_panels(input,output_dir);
  case 'layout'
   assert(ismember(j.function,{'reproduce_suppfig1','reproduce_suppfig2','reproduce_suppfig4','plot_fig2a','plot_fig2e_snr','plot_fig2e_replacement','reproduce_edfig15a','reproduce_edfig15b','plot_fig3jk_ed11b'}));
   feval(j.function,input,output_dir);
  case 'surface'
   assert(~isempty(regexp(j.stem,'^[A-Za-z][A-Za-z0-9_]*$','once')));render_brain_map(input,fullfile(output_dir,j.stem));
  otherwise,error('Unknown rendering job');
 end
end
end
