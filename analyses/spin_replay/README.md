# Replaying spatial rotations

Call `replay_spin_sources(input_file, output_file)` in MATLAB. The input MAT file contains `left_centroids` and `right_centroids` (ROI by three coordinates), and a cell array `records`. Each record contains `key`, aligned ROI vectors `x` and `y`, a MATLAB `rng` structure, and `n_rot`. ROI vectors must follow the left-then-right centroid order.

The output stores rotation indices, the observed correlation and both directional null correlations. With 1,000 rotations, the paired distribution contains 2,000 correlations. Plotting uses `mean(abs(spin_null) >= abs(observed))` without a +1 correction. These paired values are not 2,000 independently sampled rotations.
