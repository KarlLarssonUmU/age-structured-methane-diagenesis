function plotResults(r)
%PLOTRESULTS Main-model fits, source profiles, and age-class export fractions.
out=r.config.outputDir; z=r.z*100;
blue=[.08 .33 .52]; orange=[.76 .31 .12];
fig=figure('Visible','off','Color','w','Position',[100 100 1250 820]);
cleanup=onCleanup(@() close(fig));
tl=tiledlayout(fig,2,3,'TileSpacing','compact','Padding','compact');
title(tl,'Reservoir-corrected main model','FontWeight','bold');
nexttile; plot(r.concentration,z,'Color',blue,'LineWidth',1.8); hold on;
plot(r.data.concentration,r.data.concentrationDepth*100,'o','Color',orange, ...
    'MarkerFaceColor','w','MarkerSize',5);
depthAxes(); ylim([0,max(r.data.concentrationDepth)*100+5]);
xlim([0,1.05*max(r.data.concentration)]);
xlabel('CH_4 (mM)'); title('Concentration');
legend('Model','Observations','Location','southwest');
nexttile; plot(r.fm,z,'Color',blue,'LineWidth',1.8); hold on;
errorbar(r.data.fm,r.data.fmDepth*100,r.data.observedFmSd,'horizontal','o', ...
    'Color',orange,'MarkerFaceColor','w','MarkerSize',5);
depthAxes(); ylim([0,max(r.data.fmDepth)*100+5]);
xlabel('F_m'); title('Methane radiocarbon');
nexttile; plot(r.productionScaled/r.config.scaling*86400*1000,z, ...
    'Color',blue,'LineWidth',1.8); hold on;
yline(r.config.lossDepth*100,'--','Loss-layer base','Color',orange);
depthAxes(); ylim([0 80]); xlabel('Production (\mumol L_{bulk}^{-1} day^{-1})');
title(sprintf('Smooth bump: peak %.1f cm',r.productionFit.shape(1)*100));
% Duplicate the layer-base depth to display the exact discontinuity.
inside=r.z<=r.config.lossDepth;
lossZ=[z(inside);r.config.lossDepth*100;r.config.bottomDepth*100];
lossValue=[r.diagnostics.lossProfile(inside)*86400*1000;0;0];
nexttile; plot(lossValue,lossZ,'Color',orange,'LineWidth',1.8);
depthAxes(); ylim([0 20]); xlabel('Loss (\mumol L_{bulk}^{-1} day^{-1})');
title(sprintf('Linear loss: k = %.5f day^{-1}',r.kLoss_per_day));
nexttile; plot(100*r.diagnostics.cumulativeDepthFraction,z,'Color',blue,'LineWidth',1.8);
hold on; yline(10,'--','10 cm','Color',orange); depthAxes();
ylim([0,r.config.bottomDepth*100]);
xlabel('Cumulative share of dissolved export (%)'); xlim([0 100]); title('Source depth');
nexttile; plot(r.diagnostics.uniqueSourceAge,100*r.diagnostics.cumulativeAgeFraction, ...
    'Color',blue,'LineWidth',1.8); grid on; box on; ylim([0 100]);
xlabel('Source apparent radiocarbon age (yr BP)');
ylabel('Cumulative share of dissolved export (%)'); title('Source age');
exportgraphics(fig,fullfile(out,'main_model.png'),'Resolution',200);
exportgraphics(fig,fullfile(out,'main_model.pdf'),'ContentType','vector');
end

function depthAxes()
set(gca,'YDir','reverse'); ylabel('Depth (cm)'); grid on; box on;
end
