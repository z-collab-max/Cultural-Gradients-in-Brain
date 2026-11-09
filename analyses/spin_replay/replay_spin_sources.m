function replay_spin_sources(input_file,output_file)
% Recreate paired-direction spin correlations from saved rotation states.
addpath(fileparts(mfilename('fullpath')));S=load(input_file);result=struct;
for j=1:numel(S.records)
 d=S.records{j};state=d.rng;state.State=state.State(:);rng(state);
 indices=rotate_parcellation(S.left_centroids,S.right_centroids,d.n_rot);
 x=d.x(:);y=d.y(:);observed=corr(x,y);
 null=[corr(x(indices),y);corr(y(indices),x)];
 result.(d.key)=struct('spin_null',null(:),'observed',observed,'x',x,'y',y,'spin_indices_1based',indices);
end
save(output_file,'-struct','result','-v7');
end
