function run_vertex_spin(pkg,enigma,out)
if ~isfolder(out),mkdir(out);end
addpath(genpath(fullfile(enigma, 'shared')));
res = fullfile(pkg, 'vertiex_analysis_expcs', 'results');
g = load(fullfile(pkg, 'gradient_maps', 'biological_gradients_fsaverage5.mat'));
a = load(fullfile(res, 'expcs_vertex_results_thickness_pc7_1.mat'), 'coeff_val','tval_val');
b = load(fullfile(res, 'expcs_vertex_results_thickness_pc7_2.mat'), 'coeff_val','tval_val');
maps = {double(a.coeff_val(:)), double(b.coeff_val(:))};
references = {-[g.eigenmode2_lh(:); g.eigenmode2_rh(:)], ...
               [g.SAaxis_lh(:); g.SAaxis_rh(:)]};
labels = {'PC1_vs_minus_eigenmode2'; 'minus_PC2_vs_SAaxis'};
maps{2} = -maps{2};
mask_l = dlmread(fullfile(pkg, 'BrainEigenmodes', 'data', ...
    'template_surfaces_volumes', 'fsaverage5_10k_cortex-lh_mask.txt')) == 1;
mask_r = dlmread(fullfile(pkg, 'BrainEigenmodes', 'data', ...
    'template_surfaces_volumes', 'fsaverage5_10k_cortex-rh_mask.txt')) == 1;
cortex = [mask_l(:); mask_r(:)];
sl = SurfStatReadSurf1(fullfile(enigma, 'shared', 'surfaces', 'fsa5_sphere_lh'));
sr = SurfStatReadSurf1(fullfile(enigma, 'shared', 'surfaces', 'fsa5_sphere_rh'));
cl = double(sl.coord.'); cr = double(sr.coord.');
assert(size(cl,1) == 10242 && size(cr,1) == 10242);
assert(all(cellfun(@numel,maps)==20484) && all(cellfun(@numel,references)==20484));
n_rot = 1000; rng(42, 'twister');
null_r = nan(n_rot,2); observed_r = nan(2,1); n_vertices = zeros(2,1);
valid = cell(2,1);
for j = 1:2
    valid{j} = cortex & isfinite(maps{j}) & isfinite(references{j}) ...
        & maps{j} ~= 0 & references{j} ~= 0;
    observed_r(j) = corr(maps{j}(valid{j}), references{j}(valid{j}));
    n_vertices(j) = nnz(valid{j});
end
reflection = diag([-1,1,1]);
for i = 1:n_rot
    [q,rr] = qr(randn(3)); q = q * diag(sign(diag(rr)));
    if det(q) < 0, q(:,1) = -q(:,1); end
    qr_hemi = reflection * q * reflection;
    % Nearest-vertex interpolation on the sphere, separately by hemisphere.
    il = knnsearch(cl, cl*q);
    ir = knnsearch(cr, cr*qr_hemi) + 10242;
    ix = [il; ir];
    for j = 1:2
        spun = maps{j}(ix);
        ok = valid{j} & isfinite(spun) & spun ~= 0;
        if nnz(ok) > 100
            null_r(i,j) = corr(spun(ok), references{j}(ok));
        end
    end
end
valid_rotations = sum(isfinite(null_r),1).';
assert(all(valid_rotations == n_rot), 'Some rotations had too few valid vertices.');
p_spin = mean(abs(null_r) >= abs(observed_r.'),1).';
results = table(labels, observed_r, p_spin, n_vertices, valid_rotations, ...
    'VariableNames', {'comparison','pearson_r','p_spin_two_sided', ...
                      'n_vertices','n_rotations'});
writetable(results, fullfile(out, 'pairwise_vertex_spin_results.csv'));
save(fullfile(out, 'pairwise_vertex_spin_nulls.mat'), 'null_r', 'results', 'n_rot');
disp(results);

addpath(fullfile(fileparts(mfilename('fullpath')),'..','..','figures'));
for pc=1:2
    reference=references{pc};effect=maps{pc};mask=valid{pc};spin_null=null_r(:,pc);
    if pc==1,t=a.tval_val(:);else,t=-b.tval_val(:);end
    values_lh=t(1:10242);values_rh=t(10243:end);color_limits=[-5 5];
    r=observed_r(pc);p=p_spin(pc);stem=sprintf('EDFig4_PC%d',pc);file=fullfile(out,[stem '.mat']);
    save(file,'reference','effect','mask','spin_null','r','p','values_lh','values_rh','color_limits');
    render_vertex_map(file,fullfile(out,stem));
end
end
