function clean=study_china_adjust(T,Y,culture)
% Remove fitted fixed covariate effects while retaining culture and group variation.
addpath(fullfile(fileparts(mfilename('fullpath')),'..','..','china'));
names=cellstr(compose('study_culture%d',1:size(culture,2)));
for k=1:numel(names),T.(names{k})=zscore(culture(:,k));end
clean=zeros(size(Y));
for r=1:size(Y,2)
 D=T;D.response=zscore(Y(:,r));m=fitlme(D,china_formula('response',names),'FitMethod','REML');
 clean(:,r)=residuals(m,'Conditional',false);
 for k=1:numel(names)
  j=strcmp(m.Coefficients.Name,names{k});clean(:,r)=clean(:,r)+D.(names{k})*m.Coefficients.Estimate(j);
 end
end
end
