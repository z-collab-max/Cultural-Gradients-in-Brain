function render_gradients(gradients_file,roi_order_file,out)
here=fileparts(mfilename('fullpath'));addpath(fullfile(here,'..','..','figures'));
d=load(gradients_file);order=load(roi_order_file,'roi_names');roi_names=order.roi_names;
assert(size(d.maps,1)==numel(roi_names));if ~isfolder(out),mkdir(out);end
for k=1:size(d.maps,2)
    values=zscore(d.maps(:,k));color_limits=[-3 3];
    stem=sprintf('gradient_%02d',k);file=fullfile(out,[stem '.mat']);
    save(file,'roi_names','values','color_limits');render_brain_map(file,fullfile(out,stem));
end
end

