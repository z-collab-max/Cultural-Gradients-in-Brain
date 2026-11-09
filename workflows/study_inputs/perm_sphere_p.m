
function [p_perm,null_dist] = perm_sphere_p(x,y,perm_id,corr_type)

if nargin < 4
    corr_type = 'spearman';
end

nroi = size(perm_id,1);  % number of regions
nperm = size(perm_id,2); % number of permutations

rho_emp = corr(x,y,'type',corr_type);     % empirical correlation

x_perm = zeros(nroi,nperm);
y_perm = zeros(nroi,nperm);
for r = 1:nperm
    for i = 1:nroi
        x_perm(i,r) = x(perm_id(i,r)); % permute x
        y_perm(i,r) = y(perm_id(i,r)); % permute y
    end
end

rho_null_xy = zeros(nperm,1);
for r = 1:nperm
    rho_null_xy(r) = corr(x_perm(:,r),y,'type',corr_type); % correlate permuted x to unpermuted y
    rho_null_yx(r) = corr(y_perm(:,r),x,'type',corr_type); % correlate permuted y to unpermuted x
end

null_dist = [rho_null_xy(:);rho_null_yx(:)];
p_perm = mean(abs(null_dist) >= abs(rho_emp));


end

