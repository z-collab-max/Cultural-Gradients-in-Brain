function publication_style(fig)
if nargin<1,fig=gcf;end
font='Arial';font_size=8;line_width=.85;
for f=reshape(fig,1,[])
 ax=findall(f,'Type','axes');
 for a=reshape(ax,1,[])
  set(a,'FontName',font,'FontUnits','points','FontSize',font_size,'FontWeight','normal','LineWidth',line_width,'LabelFontSizeMultiplier',1);
  labels=[a.XLabel a.YLabel a.ZLabel];
  if numel(a.YAxis)>1,labels=[labels a.YAxis(1).Label a.YAxis(2).Label];end
  for h=labels
   set(h,'FontName',font,'FontUnits','points','FontSize',font_size,'FontWeight','normal');
   if ischar(h.String)||isstring(h.String)||iscell(h.String),h.String=regexprep(h.String,'\\fontsize\{[^}]*\}','');end
  end
 end
 for h=reshape(findall(f,'-property','LineWidth'),1,[])
  if ~strcmp(h.Type,'figure'),h.LineWidth=line_width;end
 end
 for h=reshape(findall(f,'Type','text'),1,[]),h.FontName=font;end
 for type={'legend','colorbar'}
  for h=reshape(findall(f,'Type',type{1}),1,[])
   set(h,'FontName',font,'FontSize',font_size,'LineWidth',line_width);
   if strcmp(type{1},'colorbar'),set(h.Label,'FontName',font,'FontSize',font_size);end
  end
 end
end
end
