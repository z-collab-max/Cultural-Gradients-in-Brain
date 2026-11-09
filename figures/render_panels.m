function render_panels(input_file,output_root)
here=fileparts(mfilename('fullpath'));root=fileparts(here);
addpath(fullfile(here,'helpers'));addpath(fullfile(root,'simulation'));
if nargin<2,output_root=fullfile(root,'outputs','figures');end
S=load(input_file,'panels');assert(isfield(S,'panels'),'Input must contain panels');
for q=1:numel(S.panels)
 d=S.panels{q};assert(~isempty(regexp(d.folder,'^[A-Za-z][A-Za-z0-9_]*$','once')));
 assert(~isempty(regexp(d.name,'^[A-Za-z][A-Za-z0-9_]*$','once')));
 fields=fieldnames(d);for k=1:numel(fields),if isnumeric(d.(fields{k})),d.(fields{k})=double(d.(fields{k}));end,end
 f=figure('Color','white','Visible','off','Position',.9*[100 100 450 250]);ax=axes(f);
 switch d.kind
 case 'combat_comparison'
  close(f);f=draw_combat_comparison(d.x(:),d.y(:),true);ax=gca;xlabel(d.xlabel,'FontSize',9,'FontWeight','bold');ylabel(d.ylabel,'FontSize',9,'FontWeight','bold');
 case 'scatter'
  if isfield(d,'stat_null'),label='partial \beta';if isfield(d,'stat_label'),label=d.stat_label;end;d.stats=sprintf('%s = %.3f, p = %.4f',label,d.stat_observed,permutation_p_absolute(d.stat_null,d.stat_observed));end
  set(f,'Position',.8*[100 200 550 300]);
  splot_hemis_bound(d.x(:),d.y(:),{'#F5CBE5','#BDE5DE'},'none',[153 112 171]/255,[153 112 171]/255,logical(d.lh(:)),[],10000,[],d.stats);
  xlabel(d.xlabel);ylabel(d.ylabel);axis square;
 case 'violin'
  b=d.beta(:);p=permutation_p_absolute(d.null,b);
  if isfield(d,'additional_null'),p=max(p,permutation_p_absolute(d.additional_null,d.additional_observed));end
  p=p(:);[~,order]=sort(abs(b),'descend');n=numel(b);
  for j=1:n,k=order(j);sgn=2*(b(k)>=0)-1;plotHorizontalViolinCorr_public(ax,sgn*d.null(:,k),sgn*b(k),p(k),n-j+1,[217 227 243]/255);end
  labels=string(d.labels(:));yticks(1:n);yticklabels(flip(labels(order)));ylim([.5 n+.5]);xlim([-1 1]);xline(0,'k--','LineWidth',.85);
  xlabel('Coefficient (\beta)','FontSize',11);box off;set(ax,'TickDir','out','LineWidth',.85,'FontSize',8,'TickLength',4*[.01 .025]);
 case 'bootstrap_histogram'
  plot_saved_bootstrap(d.null(:),d.observed,'Bootstrap',d.xlabel);
 case 'histogram'
  digits=4;if strcmp(d.folder,'Figure4'),digits=3;end
  set(f,'Position',.8*[100 200 550 300]);plot_saved_permutation(d.null(:),d.observed,'Permutation test',d.xlabel,d.p,digits);
 case 'trajectory'
  d.p=arrayfun(@(j)permutation_p_absolute(d.(sprintf('null_%d',j-1)),d.r(j)),1:numel(d.r));
  set(f,'Position',[415 496 810 284]);hold on;n=numel(d.r);
  for j=1:n
   color=[225 230 240]/255;if j==double(d.selected)+1,color=[239 191 184]/255;end
   [den,xi]=ksdensity(d.(sprintf('null_%d',j-1)));w=den/max(den)*.4;patch([j-w,j+fliplr(w)],[xi,fliplr(xi)],color,'FaceAlpha',1,'EdgeColor','none');
  end
  plot(1:n,d.r,'-o','Color',[.85 .325 .098],'LineWidth',.85,'MarkerSize',8,'MarkerFaceColor',[.85 .325 .098],'MarkerEdgeColor','none');
  for j=1:n,if d.p(j)<.05,text(j,d.r(j),'*','FontSize',15,'FontWeight','bold','HorizontalAlignment','center','VerticalAlignment','bottom');end,end
  yline(0,'k--','LineWidth',.85);xticks(1:n);xticklabels(0:n-1);if isfield(d,'components'),xticklabels(d.components);end;xlim([.5 n+.5]);ylim([-1 1]);xlabel('Number of confound PCs removed');ylabel('Spatial correlation');box off;set(ax,'TickDir','out','FontSize',8,'LineWidth',.85);
 case 'vertical_violin'
  d.modes=double(d.modes);[~,order]=sort(abs(d.beta),'descend');
  hold on;for j=1:numel(order),k=order(j);sgn=2*(d.beta(k)>=0)-1;[den,xi]=ksdensity(sgn*d.null_beta(:,k));w=den/max(den)*.35;patch([j-w,j+fliplr(w)],[xi,fliplr(xi)],[217 227 243]/255,'EdgeColor','none');plot(j,abs(d.beta(k)),'o','MarkerFaceColor',[113 129 158]/255,'MarkerEdgeColor','none','MarkerSize',7);end
  xticks(1:numel(order));xticklabels(string(d.modes(order)));xlabel('Eigenmode');ylabel('Coefficient (\beta)');box off;set(ax,'TickDir','out','FontSize',8,'LineWidth',.85);
 case 'pca_biplot'
  set(f,'Position',[100 100 520 420]);hold on;scatter(d.scores(:,1),d.scores(:,2),12,[113 129 158]/255,'filled','MarkerEdgeColor','none');
  scale=.7*max(abs(d.scores),[],'all')/max(abs(d.loadings),[],'all');
  for j=1:size(d.loadings,1),v=d.loadings(j,:)*scale;plot([0 v(1)],[0 v(2)],'Color',[184 128 118]/255,'LineWidth',.85);text(v(1)*1.1,v(2)*1.1,string(d.labels(j)),'HorizontalAlignment','center','FontSize',8,'Color',[171 77 62]/255);end
  xline(0,'--','Color',[.4 .4 .4]);yline(0,'--','Color',[.4 .4 .4]);xlabel(sprintf('PC1 (%.2f%%)',d.explained(1)));ylabel(sprintf('PC2 (%.2f%%)',d.explained(2)));axis equal;box off;
 case 'identity_scatter'
  set(f,'Position',[366 450 320 280]);x=d.x(:);y=d.y(:);lim=1.12*max(abs([x;y]));lims=[-lim lim];
  scatter(x,y,12,[118 42 131]/255,'filled');hold on;plot(lims,lims,'r--','LineWidth',1);xlim(lims);ylim(lims);axis square;
  if isfield(d,'fit_line')&&d.fit_line,xx=[min(x),max(x)];plot(xx,polyval(polyfit(x,y,1),xx),'k-','LineWidth',.85);end
  stats=sprintf('r = %.3f',corr(x,y));if isfield(d,'stats'),stats=d.stats;end
  text(lims(1)+.06*range(lims),lims(2)-.12*range(lims),stats,'FontSize',9,'FontWeight','bold');xlabel(d.xlabel,'FontSize',9,'FontWeight','bold');ylabel(d.ylabel,'FontSize',9,'FontWeight','bold');box off;set(ax,'TickDir','out','LineWidth',.85,'FontSize',8,'TickLength',3*[.01 .025]);
 case 'loo'
  set(f,'Position',[100 200 650 250]);hold on;v=d.y(:);x=1:numel(v);col=[85 122 165]/255;if contains(d.name,'PC2'),col=[78 146 126]/255;end
  for k=1:numel(v),plot([x(k) x(k)],[v(k) 1],'-','Color',.22*col+.78*[1 1 1],'LineWidth',1.15);end
  scatter(x,v,42,col,'filled','MarkerEdgeColor','w','LineWidth',.6);yline(1,'--','Color',[.6 .6 .6],'LineWidth',.8);
  set(ax,'XLim',[.3 numel(v)+.7],'YLim',[min(.55,min(v)-.03) 1.015],'YTick',.6:.1:1,'XTick',x,'XTickLabel',string(d.labels),'XTickLabelRotation',58,'TickLabelInterpreter','none','FontName','Arial','FontSize',9,'Box','off','YGrid','on','GridAlpha',.12);ylabel('Spatial correlation with full map','FontSize',11);
 case 'curves'
  hold on;colors=[184 128 118;113 129 158]/255;
  if isfield(d,'colors'),colors=d.colors;end
  for j=1:size(d.y,2),ci=min(j,size(colors,1));plot(d.x,d.y(:,j),'-o','Color',colors(ci,:),'LineWidth',1.5,'MarkerSize',4,'MarkerFaceColor',colors(ci,:));end
  legend(string(d.labels),'Location','best','Box','off');xlabel(d.xlabel);ylabel(d.ylabel);box off;set(ax,'TickDir','out','FontSize',8,'LineWidth',.85);
 case 'ribbons'
  hold on;colors=[.85 .33 .10;0 .45 .74];
  for j=1:size(d.series,1)
   a=squeeze(d.series(j,:,:));xx=d.x(:);ok=all(isfinite(a),2);xx=xx(ok);a=a(ok,:);c=colors(min(j,2),:);
   if numel(d.labels)==1&&strcmp(string(d.labels),'Subject'),c=colors(2,:);end
   fill([xx;flipud(xx)],[a(:,2);flipud(a(:,3))],c,'FaceAlpha',.2,'EdgeColor','none','HandleVisibility','off');plot(xx,a(:,1),'-o','Color',c,'LineWidth',1.5,'MarkerSize',4,'MarkerFaceColor',c);
  end
  legend(string(d.labels),'Location','best','Box','off');xlabel('Sampling proportion');ylabel('Spatial correlation');box off;set(ax,'TickDir','out','FontSize',8,'LineWidth',.85);
 case 'heatmap'
  if isfield(d,'cell_nulls'),d.p=zeros(size(d.matrix));for ii=1:size(d.matrix,1),for jj=1:size(d.matrix,2),d.p(ii,jj)=permutation_p_absolute(d.cell_nulls(:,ii,jj),d.matrix(ii,jj));end,end,end
  if isfield(d,'spin_null')&&isfield(d,'permutation_null'),d.p=heatmap_p_absolute(d);end
  set(f,'Position',[100 100 600 520]);
  if isfield(d,'zero_diagonal') && d.zero_diagonal && size(d.matrix,1)==size(d.matrix,2),d.matrix(1:size(d.matrix,1)+1:end)=0;end
  imagesc(d.matrix,d.clim);axis image;
  colormp_1;cm=mycolor;if isfield(d,'colors'),cm=d.colors;end;colormap(ax,cm);colorbar;
  if isfield(d,'split'),hold on;xline(d.split,'k-');yline(d.split,'k-');end
  labs=string(d.labels);xl=labs;if isfield(d,'xlabels'),xl=string(d.xlabels);end
  xticks(1:size(d.matrix,2));xticklabels(xl);yticks(1:size(d.matrix,1));yticklabels(labs);xtickangle(45);set(ax,'FontSize',8,'TickLength',[0 0]);
  if ~startsWith(d.name,'Fig4_')
   for ii=1:size(d.matrix,1),for jj=1:size(d.matrix,2)
    txt=sprintf('%.2f',d.matrix(ii,jj));stars='';if isfield(d,'p'),pv=d.p(ii,jj);if pv<.001,stars='***';elseif pv<.01,stars='**';elseif pv<.05,stars='*';end,end
    tc='k';if abs(d.matrix(ii,jj))>.6,tc='w';end
    if size(d.matrix,1)~=size(d.matrix,2)||ii~=jj,text(jj,ii,txt,'Color',tc,'FontSize',6,'HorizontalAlignment','center');text(jj,ii-.25,stars,'Color',tc,'FontSize',7,'FontWeight','bold','HorizontalAlignment','center');end
   end,end
  end
 case 'bars'
  if ~isfield(d,'color'),d.color=repmat([.44 .51 .62],numel(d.y),1);end
  set(f,'Position',[100 100 650 260]);b=bar(d.y,.85,'FaceColor','flat','EdgeColor','none');b.CData=d.color;xticks(1:numel(d.y));xticklabels(string(d.labels));if numel(d.y)>10,xtickangle(90);end;ylabel(d.ylabel);box off;set(ax,'TickDir','out','LineWidth',.85,'FontSize',8);
 case 'error_scatter'
  set(f,'Position',[100 100 380 320]);x=d.x(:);y=d.y(:);se=d.se(:);hold on;
  if strcmp(d.folder,'ExtendedDataFigure1')
   errorbar(x,y,se,se,d.xse(:),d.xse(:),'Color','#E18F81','LineStyle','none','LineWidth',.75);plot(x,y,'.','Color','#E18F81','MarkerSize',20);text(x,y,string(d.labels),'Color','#E18F81','FontSize',5);
   if endsWith(d.name,'c'),co=polyfit(x,y,1);xx=xlim;plot(xx,polyval(co,xx),'Color','#E18F81','LineWidth',1.5);end
   xlabel(d.xlabel);ylabel(d.ylabel);box off;set(ax,'FontName','Helvetica','FontSize',8,'TickDir','out','LineWidth',.85);
  else
  if isfield(d,'xse'),errorbar(x,y,se,se,d.xse(:),d.xse(:),'o','Color',d.color,'MarkerSize',3,'LineWidth',.8,'CapSize',2);else,errorbar(x,y,se,'o','Color',d.color,'MarkerSize',3,'LineWidth',.7,'CapSize',0);end
  coef=polyfit(x,y,1);xx=linspace(min(x),max(x),100)';A=[x ones(size(x))];covb=inv(A'*A)*sum((y-A*coef(:)).^2)/(numel(y)-2);V=[xx ones(size(xx))];yy=V*coef(:);delta=1.96*sqrt(sum((V*covb).*V,2));plot(xx,yy,'Color',[.5 .5 .5],'LineWidth',1);fill([xx;flipud(xx)],[yy-delta;flipud(yy+delta)],[.5 .5 .5],'FaceAlpha',.15,'EdgeColor','none');
  if isfield(d,'labels'),text(x,y,string(d.labels),'FontSize',5,'Color',d.color);end
  text(.04,.04,d.stats,'Units','normalized','FontSize',8);xlabel(d.xlabel);ylabel(d.ylabel);box off;set(ax,'FontSize',9,'TickDir','out','LineWidth',.85);
  end
 case 'purple_fit'
  set(f,'Position',[100 100 340 240]);hold on;x=d.x(:);y=d.y(:);mn=min([x;y]);mx=max([x;y]);rr=mx-mn;lims=[mn-.1*rr,mx+.1*rr];
  scatter(x,y,12,[128 115 172]/255,'filled','MarkerFaceAlpha',1);plot(lims,lims,'--','Color',[.6 .6 .6],'LineWidth',.85);
  coef=polyfit(x,y,1);xx=linspace(lims(1),lims(2),100);plot(xx,polyval(coef,xx),'-','Color',[.2 .2 .2],'LineWidth',.85);xlim(lims);ylim(lims);
  xlabel(d.xlabel);ylabel(d.ylabel);r=corr(x,y);tc='k';if abs(r)>.9,tc=[0 .5 0];end
  text(lims(1)+.05*range(lims),lims(2)-.15*range(lims),sprintf('r = %.3f',r),'Color',tc,'FontSize',10,'FontWeight','bold','VerticalAlignment','top');grid off;box off;set(ax,'TickDir','out','LineWidth',.85,'FontSize',8,'TickLength',4*[.01 .025]);axis square;
  if isfield(d,'p'),if d.p<.001,pt='p < 0.001';else,pt=sprintf('p = %.3f',d.p);end;text(lims(1)+.05*range(lims),lims(2)-.25*range(lims),pt,'Color',tc,'FontSize',10,'FontWeight','bold');end
 case 'task_scatter'
  set(f,'Position',[100 100 360 280]);x=double(d.x(:));y=d.y(:);hold on;mdl=fitlm(x,y);xx=linspace(min(x),max(x),100)';[yy,ci]=predict(mdl,xx);fill([xx;flipud(xx)],[ci(:,1);flipud(ci(:,2))],[.9 .9 .9],'EdgeColor','none','FaceAlpha',.5);plot(xx,yy,'Color',[.4 .4 .4],'LineWidth',.85);scatter(x,y,25,d.color,'filled');xlim([min(x)-.5 max(x)+.5]);xticks(x);xticklabels(string(d.labels));xtickangle(45);ylabel('Task effect');box off;set(ax,'FontSize',8,'TickDir','out','LineWidth',.85);
  if isfield(d,'stats'),text(.04,.94,d.stats,'Units','normalized','FontSize',8,'VerticalAlignment','top');end
 case 'paired_histogram'
  set(f,'Position',[100 100 550 300]);hold on;
  histogram(d.null_H,55,'Normalization','probability','FaceColor',[230 160 148]/255,'FaceAlpha',.75,'EdgeColor','none');histogram(d.null_L,55,'Normalization','probability','FaceColor',[141 161 198]/255,'FaceAlpha',.75,'EdgeColor','none');
  colors=[181 79 72;65 99 160;230 160 148;141 161 198]/255;values=d.means([5 6 3 4]);for j=1:4,xline(values(j),'Color',colors(j,:),'LineWidth',1.5);end
  xlim([-.2 .2]);xlabel('Positive ROIs - Negative ROIs');ylabel('Proportion');box off;set(ax,'FontSize',8,'TickDir','out','LineWidth',.85);
  inset=axes(f,'Position',[.74 .61 .20 .26]);histogram(d.null_difference,40,'Normalization','pdf','FaceColor',[.8 .8 .8],'EdgeColor','none');hold on;xline(d.cross,'Color',[.698 .094 .169]);xline(d.same,'k--');title(sprintf('p = %.3f',d.p_cross),'FontSize',8);set(inset,'FontSize',6,'YTick',[]);box off;
 case 'roi_intervals'
  set(f,'Position',[100 100 280 450]);y=(numel(d.stable):-1:1)';hold on;m=d.stable(:);v=d.mover(:);fill([m;flipud(v)],[y;flipud(y)],d.light,'FaceAlpha',.25,'EdgeColor','none');errorbar(m,y,m-d.lo(:),d.hi(:)-m,'horizontal','o','Color',d.color,'MarkerSize',3,'CapSize',2);plot(v,y,'o-','Color',d.light,'MarkerSize',3,'LineWidth',.7);xline(0,'--','Color',[.6 .6 .6]);xlim([-.2 .2]);xlabel('Z-scored thickness');set(ax,'YTick',[],'FontSize',8,'TickDir','out','LineWidth',.85);box off;
 case 'intervals'
  hold on;y=(1:numel(d.estimate))';
  errorbar(d.estimate(:),y,d.estimate(:)-d.low(:),d.high(:)-d.estimate(:),'horizontal','o','Color',[.35 .45 .65],'CapSize',3);
  yticks(y);yticklabels(string(d.labels));ylim([.5 numel(y)+.5]);xline(0,'--');xlabel(d.xlabel);box off;
 case 'surface'
  patch(ax,'Vertices',double(d.vertices),'Faces',double(d.faces),'FaceVertexCData',double(d.values(:)),'FaceColor','interp','EdgeColor','none');
  axis equal off;view(d.view);colormap(ax,parula);clim(d.clim);colorbar;
 case 'table'
  axis off;rows=string(d.labels);cols=string(d.columns);a=d.matrix;
  for j=1:numel(cols),text(j,0,cols(j),'Interpreter','none','HorizontalAlignment','center');end
  for i=1:numel(rows),text(0,-i,rows(i),'Interpreter','none');for j=1:numel(cols),text(j,-i,sprintf('%.4g',a(i,j)),'HorizontalAlignment','center');end,end
  xlim([-.5 numel(cols)+.5]);ylim([-numel(rows)-.5 .5]);
  otherwise,error('Unsupported panel kind: %s',d.kind);
 end
 out=fullfile(output_root,d.folder);if ~isfolder(out),mkdir(out);end
 publication_style(f);savefig(f,fullfile(out,[d.name '.fig']));
 exportgraphics(f,fullfile(out,[d.name '.png']),'Resolution',150);close(f);
 fprintf('EXPORTED %s/%s\n',d.folder,d.name);
end
end

