function [T,metadata]=prepare_macro_pcs(T,output_file)
% Fit once per geographic unit and join scores back in the original row order.
city={'absolute_latitude','terrain_ruggedness','historical_humidity','historical_temperature',...
 'historical_sunshine','historical_precipitation','population_density','green_space','natural_disaster_frequency_per_capita'};
province={'hdi_education','hdi_income'};
[T,metadata.city]=block(T,{'province_birth','city_birth'},city,'city_macro');
[T,metadata.province]=block(T,{'province_birth'},province,'province_macro');
if nargin>1
 parent=fileparts(output_file);if ~isempty(parent)&&~isfolder(parent),mkdir(parent);end
 save(output_file,'T','metadata');
end
end
function [T,m]=block(T,keys,fields,prefix)
assert(~any(startsWith(T.Properties.VariableNames,[prefix '_pc'])),'Use raw inputs without existing macro PCs');
assert(all(ismember([keys fields],T.Properties.VariableNames)),'Missing raw macro inputs');
assert(~any(ismissing(T(:,keys)),'all'),'Missing geographic keys');
[g,units]=findgroups(T(:,keys));X=double(T{:,fields});assert(all(isfinite(X),'all'));
assert(height(units)>=2,'PCA requires at least two geographic units');
U=zeros(height(units),numel(fields));
for j=1:height(units)
 v=X(g==j,:);assert(all(v==v(1,:),'all'),'Macro inputs vary within a geographic unit');U(j,:)=v(1,:);
end
mu=mean(U,1);scale=std(U,0,1);scale(scale==0)=1;Z=(U-mu)./scale;
[~,S,V]=svd(Z,'econ');s=diag(S);tolerance=max(size(Z))*eps*s(1);k=sum(s>tolerance);V=V(:,1:k);
for j=1:k,[~,i]=max(abs(V(:,j)));if V(i,j)<0,V(:,j)=-V(:,j);end,end
scores=Z*V;names=arrayfun(@(j)sprintf('%s_pc%d',prefix,j),1:k,'UniformOutput',false);
T(:,fields)=[];for j=1:k,T.(names{j})=scores(g,j);end
m=struct('fields',{fields},'keys',{keys},'units',units,'mean',mu,'scale',scale,'loadings',V,...
 'rank',k,'names',{names},'tolerance',tolerance,'explained',s(1:k).^2/sum(s.^2),'scores',scores);
end
