function cfg = defaultConfig()
%DEFAULTCONFIG Main manuscript model. All six free parameters are fitted.
cfg.repoDir = fileparts(fileparts(mfilename('fullpath')));
cfg.ageScenario = 'reservoirCorrected';
cfg.fmTarget = 'observations';
cfg.productionFamily = 'bump';
cfg.bottomDepth = 3.62;                 % L [m]
cfg.meshStep = 0.0025;                  % [m]; includes the loss boundary
cfg.scaling = 1e8;                     % numerical scaling only
cfg.D0 = 1.0158e-9;                     % molecular diffusion [m^2/s]
cfg.concentrationFitPoints = 8;         % measured near-surface regression
cfg.boundaryMode = 'gradient';
cfg.fixedKsw = NaN;                    % used only by explicit sensitivity runs
cfg.fitDarcy = true;
cfg.qMagnitude = 0;                    % ignored when fitDarcy=true
cfg.qBounds = [0,1e-8];                 % positive UPWARD Darcy flux [m/s]
cfg.qGrid = [0,logspace(-11,-8,17)];
cfg.fitLoss = true;
cfg.kLoss = 0;                         % ignored when fitLoss=true
cfg.kLossBounds = [0,0.25];             % [day^-1], coefficient of BULK loss
cfg.kLossGrid = [0,.001,.003,.006,.012,.025,.05,.1,.15,.25];
cfg.lossDepth = 0.10;                  % loss layer [0,z_loss], [m]
cfg.deepLossMultiplier = 0;
% Columns: peak depth and width of production tail [m]. The effective
% peak lower bound is always max(bumpBounds(1,1),lossDepth).
cfg.bumpBounds = [0,0.10;0.60,0.60];
cfg.shapeStarts = 3;
cfg.shapeTolerance = 1e-8;
cfg.outerTolerance = 1e-5;
cfg.outerStarts = 3;
cfg.outerMaxIterations = 400;
cfg.outerObjective = 'Fm';
cfg.verbose = true;
cfg.outputDir = fullfile(cfg.repoDir,'results','main');
end
