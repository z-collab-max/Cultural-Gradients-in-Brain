function result=refit_sample(T,cortex,rows)
% rows indexes the original individual table and cortical matrix.
assert(numel(rows)==numel(unique(rows)),'Use distinct rows for a deletion or subsample fit.');
D=prepare_china_table(T(rows,:));result=fit_china_mixed(D,cortex(rows,:));
end
