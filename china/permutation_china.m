function result=permutation_china(T,statistic,repetitions,seed)
% statistic accepts an input table and returns a finite numeric vector.
if nargin<3,repetitions=10000;end
if nargin<4,seed=20261103;end
assert(repetitions>0&&repetitions==floor(repetitions));
stream=RandStream('mt19937ar','Seed',seed);G=numel(unique(T.province_birth));
orders=zeros(repetitions,G);for b=1:repetitions,orders(b,:)=randperm(stream,G);end
observed=statistic(T);assert(all(isfinite(observed),'all'));observed=observed(:)';
null=zeros(repetitions,numel(observed));
for b=1:repetitions
 shuffled=shuffle_china_culture(T,orders(b,:));
 value=statistic(shuffled);null(b,:)=value(:)';
 assert(all(isfinite(null(b,:))));
end
result=struct('observed',observed,'null',null,'orders',orders,...
 'p_two_sided_absolute',mean(abs(null)>=abs(observed),1));
end
