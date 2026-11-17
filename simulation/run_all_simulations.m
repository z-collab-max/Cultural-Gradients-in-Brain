function run_all_simulations(smoke)
if nargin<1,smoke=false;end
here=fileparts(mfilename('fullpath'));old=getenv('V3_SMOKE');
cleanup=onCleanup(@()setenv('V3_SMOKE',old));
% Artifact-score rank is bounded by the number of groups.
pca_warning=warning('off','stats:pca:ColRankDefX');
restore_warning=onCleanup(@()warning(pca_warning));
if smoke,setenv('V3_SMOKE','1');else,setenv('V3_SMOKE','0');end
evalin('base',sprintf('run(''%s'')',fullfile(here,'generate_ed2.m')));
evalin('base',sprintf('run(''%s'')',fullfile(here,'generate_si3.m')));
end
