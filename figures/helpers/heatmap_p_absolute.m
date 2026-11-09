function p=heatmap_p_absolute(d)
% Null arrays have dimensions draws x rows x columns.
n=size(d.matrix,1);p=nan(n);ps=nan(n);pp=nan(n);
for i=1:n,for j=1:n
 a=d.spin_null(:,i,j);b=d.permutation_null(:,i,j);
 if all(isfinite(a)),ps(i,j)=permutation_p_absolute(a,d.matrix(i,j));end
 if all(isfinite(b)),pp(i,j)=permutation_p_absolute(b,d.matrix(i,j));end
end,end
if d.symmetric_spin,ps=max(ps,ps');end
if d.symmetric_permutation,pp=max(pp,pp');end
p=max(ps,pp);missing=isnan(pp);
if d.spin_only_cells,p(missing)=ps(missing);else,p(missing)=NaN;end
p(1:n+1:end)=NaN;
end
