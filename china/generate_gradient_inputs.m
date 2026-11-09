function generate_gradient_inputs(input_root,output_file)
% Assemble ROI gradients from hemisphere summary files and surface eigenmodes.
L=load(fullfile(input_root,'DK308_summary_L.mat'),'Summary_DK308_L');
R=load(fullfile(input_root,'DK308_summary_R.mat'),'Summary_DK308_R');
L=L.Summary_DK308_L;R=R.Summary_DK308_R;
fields={'SAaxis','intervar','cogPC1','xu_evo','reardon_dev','xu_FChomology','gene_PC1'};
S=struct;
for j=1:numel(fields)
    k=fields{j};S.(k)=[L.(k)(2:end,1);R.(k)(2:end,1)];
end
S.myelin=[L.T1w_T2w_T1wT2w(2:end,1);R.T1w_T2w_T1wT2w(2:end,1)];
S.fcgradientsM=[L.fcgradient(2:end,:);R.fcgradient(2:end,:)];
A=load(fullfile(input_root,'surface_500_L_164k_midthickness_summary.mat'),'surface_500_L');
B=load(fullfile(input_root,'surface_500_R_164k_midthickness_summary.mat'),'surface_500_R');
modes=cell(1,2);
for h=1:2
    if h==1,s=A.surface_500_L;else,s=B.surface_500_R;end
    labels=unique(s.label);m=nan(numel(labels),size(s.modes,2));
    for r=1:numel(labels)
        for k=1:size(m,2)
            idx=s.label==labels(r) & ~isnan(s.modes(:,k));
            if h==1,idx=idx & s.cortex~=0;else,idx=idx & s.label~=0;end
            m(r,k)=mean(s.modes(idx,k));
        end
    end
    modes{h}=m(2:end,:);
end
S.DK308_modes_asym=[-modes{1};modes{2}];S.mode2=S.DK308_modes_asym(:,2);
save(output_file,'-struct','S','-v7');
end
