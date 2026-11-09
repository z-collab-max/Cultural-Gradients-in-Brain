function p=permutation_p_absolute(null_values,observed)
% Absolute two-sided empirical tails, inclusive ties, without a +1 correction.
observed=double(observed(:)');null_values=double(null_values);
if isscalar(observed),null_values=null_values(:);end
assert(ismatrix(null_values)&&size(null_values,2)==numel(observed),'Null must be draws by statistics');
assert(size(null_values,1)>0&&all(isfinite(null_values),'all')&&all(isfinite(observed)),'Finite draws and observations are required');
p=mean(abs(null_values)>=abs(observed),1);
end
