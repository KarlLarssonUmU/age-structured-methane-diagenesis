function r = reactionSettings(cfg)
%REACTIONSETTINGS Linear loss: R_loss = k_loss f(z) C, in bulk-volume units.
r.fit=cfg.fitLoss; r.value=cfg.kLoss; r.bounds=cfg.kLossBounds;
r.grid=cfg.kLossGrid; r.scale=.01; r.depth=cfg.lossDepth;
r.deepMultiplier=cfg.deepLossMultiplier; r.perDayFactor=1;
r.historyName='kLoss_per_day';
end
