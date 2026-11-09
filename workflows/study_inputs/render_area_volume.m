function render_area_volume(source,out)
addpath(fullfile(fileparts(mfilename('fullpath')),'..','..','figures'));
for modality=["area","volume"]
 for pc=1:2
  stem=modality+"_PC"+pc;
  render_panels(fullfile(source,stem,'ExtendedDataFigure7','panels.mat'),out);
  render_brain_map(fullfile(source,stem+"_surface.mat"),char(fullfile(out,'ExtendedDataFigure7',stem+"_surface")));
 end
end
end
