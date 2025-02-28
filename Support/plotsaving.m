FolderName = '.plots';   % Your destination folder
FigList = findobj(allchild(0), 'flat', 'Type', 'figure');
for iFig = 1:length(FigList)
  FigHandle = FigList(iFig);
  %fontsize(FigHandle, 11, "points") % fontsize
  FigName   = num2str(get(FigHandle, 'Number'));
  set(0, 'CurrentFigure', FigHandle);
  %--- reorder stacking of lines
  %set(gca, 'Children', flipud(get(gca, 'Children')) ) 
  savefig(fullfile(FolderName, [FigName '.fig']));
  %--- save with black background
  %exportgraphics(gcf, fullfile(FolderName, [FigName '.pdf']),'BackgroundColor','k','Resolution',600)
  %exportgraphics(gcf, fullfile(FolderName, [FigName '.png']),'BackgroundColor','k','Resolution',600)
  %--- save with white background
  exportgraphics(gcf, fullfile(FolderName, [FigName '.pdf']),'Resolution',600)
  exportgraphics(gcf, fullfile(FolderName, [FigName '.png']),'Resolution',600)
end