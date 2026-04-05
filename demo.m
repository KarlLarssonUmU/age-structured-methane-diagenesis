function demo()
%DEMO Run a minimal age-structured methane diagenetic model example.
%   This script reproduces the example figures and spreadsheet in the
%   results/ folder for the public repository. The Darcy flux q is fixed,
%   while the methane production profile parameters are optimized against
%   porewater CH4 concentration and Fm observations.

% Make helper functions available regardless of the current working folder.
repoDir = fileparts(mfilename('fullpath'));
addpath(fullfile(repoDir, 'src'));

%% User-configurable settings
qDarcy = -8.198032e-10;     % Darcy flux [m/s]; negative means upward flow
bottomDepth = 3.62;         % Lower boundary depth [m]
betaF = 0;                  % Weight on the Fm misfit in the inverse problem
scaling = 1e8;              % Numerical scaling used to improve conditioning
concentrationFitPoints = 8; % Number of shallow data points used to fit k_sw

%% Load input data
dataFile = fullfile(repoDir, 'data', 'data.xlsx');
data = readtable(dataFile, 'Sheet', 1);

% Sediment porosity profile
porosityDepth = data.poro_cm / 100;    % [m]
porosity = data.poro_perc / 100;       % fraction [-]
idx = ~isnan(porosity);
porosityDepth = porosityDepth(idx);
porosity = porosity(idx);

% Sediment radiocarbon signature profile
sedimentDepth = data.sedi_cm / 100;    % [m]
sedimentFm = data.sedi_Fm;             % Fm [-]
sedimentAgeBP = data.sedi_age;         % years BP
idx = ~isnan(sedimentFm);
sedimentDepth = sedimentDepth(idx);
sedimentFm = sedimentFm(idx);
sedimentAgeBP = sedimentAgeBP(idx);

% Methane Fm observations
methaneFmDepth = data.ch4a_cm / 100;   % [m]
methaneFm = data.ch4a_Fm;              % Fm [-]
idx = ~isnan(methaneFm);
methaneFmDepth = methaneFmDepth(idx);
methaneFm = methaneFm(idx);

% Methane concentration observations
methaneConcDepth = data.ch4c_cm / 100; % [m]
methaneConc = data.ch4c_conc;          % [mM]
idx = ~isnan(methaneConc);
methaneConcDepth = methaneConcDepth(idx);
methaneConc = methaneConc(idx);

%% Transport parameters
% Tortuosity following Boudreau (1996).
theta2 = 1 - log(porosity.^2);

% Molecular diffusion coefficient of CH4 in free water [m^2/s].
D0 = 1.0158e-9;

% Effective molecular diffusion coefficient in porewater [m^2/s].
D_eff = D0 ./ theta2;

% Full diffusion coefficient phi * D_eff, scaled for numerical stability.
coeffData = scaling * D_eff .* porosity;
coeffFcn = @(depth) interp1(porosityDepth, coeffData, depth, 'linear', 'extrap');

% Fit a straight line to the shallowest concentration points to estimate
% the sediment-water transfer coefficient k_sw.
lineFit = polyfit(methaneConcDepth(1:concentrationFitPoints), methaneConc(1:concentrationFitPoints), 1);
dCdz0 = lineFit(1);
C0 = lineFit(2);

%% Finite-element discretization
topDepth = 0;
nElements = round(100 * (bottomDepth - topDepth)); % 1 cm resolution
basisType = 'p1';
reactionFcn = @(z) 0; % no methane consumption in this minimal example

%% Assemble transport operator for the prescribed Darcy flux
qDarcyScaledFcn = @(x) qDarcy * scaling;
kSwScaled = coeffFcn(0) * dCdz0 / C0 - qDarcyScaledFcn(0);

[A, M, ~, meshNodes, meshElements, basisCoords, basisConn] = assemble1dAdvectionDiffusion( ...
    topDepth, bottomDepth, nElements, basisType, coeffFcn, reactionFcn, @(x) 0, qDarcyScaledFcn);
[A_bc, ~] = assemble1dBC(meshNodes, meshElements, basisCoords, basisConn, basisType, kSwScaled, 0);

% Matrix mapping a production source vector to the solution vector.
responseMatrix = (A + A_bc) \ M;

%% Solve the inverse problem for the production profile
[productionCoeffs, E_tot, E_C, E_F, productionMap, ~, ~, ~] = optimizeProductionProfile( ...
    basisCoords, meshNodes, responseMatrix, sedimentDepth, sedimentFm, ...
    methaneConcDepth, methaneConc, methaneFmDepth, methaneFm, betaF);

% Modeled concentration profile at all mesh nodes.
modeledConcentration = responseMatrix * productionMap' * productionCoeffs;

% Modeled production profile at all mesh nodes (scaled units).
modeledProductionScaled = productionMap' * productionCoeffs;

% In this minimal example, methane produced at a given depth is assumed to
% have the sediment Fm corresponding to that same depth.
sedimentFmFcn = @(z) interp1(sedimentDepth, sedimentFm, z, 'linear', 'extrap');
Fm_nodes = sedimentFmFcn(meshNodes);
modeledMeanFm = (responseMatrix * diag(Fm_nodes) * modeledProductionScaled) ./ ...
                (responseMatrix * modeledProductionScaled);

fprintf('q = %.2e m/s: E_tot=%.3e (E_C=%.3e, E_F=%.3e)\n', qDarcy, E_tot, E_C, E_F);

%% Ensure the output folder exists
resultsDir = fullfile(repoDir, 'results');
if ~exist(resultsDir, 'dir')
    mkdir(resultsDir);
end

%% Plot concentration profile
figWidth = 200; figHeight = 470;
labels = {'Simulation', 'Data'};

fig = figure(1); clf; hold on
fig.Units = 'pixels';
fig.Position(3:4) = [figWidth, figHeight];
plot(modeledConcentration, meshNodes, 'LineWidth', 3);
plot(methaneConc, methaneConcDepth, 'kx', 'LineWidth', 1.0);
set(gca, 'YDir', 'reverse')
xlabel('Concentration (mM)')
ylabel('Depth (m)')
legend(labels, 'Location', 'south')
exportgraphics(gca, fullfile(resultsDir, 'concentration-profile.png'));

%% Plot mean Fm profile
fig = figure(2); clf; hold on
fig.Units = 'pixels';
fig.Position(3:4) = [figWidth, figHeight];
plot(modeledMeanFm, meshNodes, 'LineWidth', 3);
plot(methaneFm, methaneFmDepth, 'kx', 'LineWidth', 0.5);
set(gca, 'YDir', 'reverse')
xlabel('Radiocarbon signature (F_m)')
ylabel('Depth (m)')
legend(labels, 'Location', 'south')
exportgraphics(gca, fullfile(resultsDir, 'mean-radiocarbon-age.png'));

%% Plot methane production profile
fig = figure(3); clf; hold on
fig.Units = 'pixels';
fig.Position(3:4) = [figWidth, figHeight];
plot(modeledProductionScaled, meshNodes, 'LineWidth', 3);
set(gca, 'YDir', 'reverse')
xlabel('Production [10^{-8} mol m^{-3} s^{-1}]')
ylabel('Depth [m]')
legend({'Simulation'}, 'Location', 'south')
if exist('ysecondarylabel', 'file') == 2
    ysecondarylabel(Visible="off")
end
exportgraphics(gca, fullfile(resultsDir, 'production-shape.png'));

%% Compute cumulative age distribution of the upward CH4 flux
topSampleMatrix = sampleMatrixP1(0, meshNodes);
ageSpecificProductionScaled = productionMap' * productionCoeffs;
ageSpecificConcentrationAtTop = zeros(size(ageSpecificProductionScaled));
for i = 1:length(ageSpecificProductionScaled)
    ageSpecificConcentrationAtTop(i) = topSampleMatrix * responseMatrix(:, i) * ageSpecificProductionScaled(i);
end

totalTopConcentration = sum(ageSpecificConcentrationAtTop);
sedimentAgeClassesBP = round(sedimentAgeBP(1:length(ageSpecificProductionScaled)));
ageSpecificConcentrationAtTop(end) = 2 * ageSpecificConcentrationAtTop(end);

relativeFlux = ageSpecificConcentrationAtTop / totalTopConcentration;
cumulativeRelativeFlux = cumsum(relativeFlux);

ageGridBP = sedimentAgeClassesBP;
cumulativeFluxFcn = @(age) interp1(ageGridBP, cumulativeRelativeFlux, age, 'linear', 'extrap');

targetFractions = [0.10, 0.25, 0.50, 0.75, 0.90];
ageQuery = linspace(min(ageGridBP), max(ageGridBP), 1000);
cumulativeValues = arrayfun(cumulativeFluxFcn, ageQuery);
targetAgesBP = zeros(size(targetFractions));
for k = 1:numel(targetFractions)
    tf = targetFractions(k);
    if tf <= min(cumulativeValues)
        targetAgesBP(k) = ageQuery(1);
    elseif tf >= max(cumulativeValues)
        targetAgesBP(k) = ageQuery(end);
    else
        targetAgesBP(k) = interp1(cumulativeValues, ageQuery, tf, 'linear');
    end
end

positiveAgeMask = find(ageGridBP >= 1);
xVals = [1; ageGridBP(positiveAgeMask)];
yVals = [cumulativeFluxFcn(1); 100 * cumulativeRelativeFlux(positiveAgeMask)];

fig = figure(4); clf;
semilogx(xVals, yVals, 'LineWidth', 3)
fluxLegend = {'Cumulative %'};
ylim([0 100]);
hold on
yLimits = ylim;
for k = 1:numel(targetAgesBP)
    x = targetAgesBP(k);
    if k == 3
        semilogx([x x], yLimits, 'r-', 'LineWidth', 1.5);
    elseif k == 2 || k == 4
        semilogx([x x], yLimits, 'k-', 'LineWidth', 0.75);
    else
        semilogx([x x], yLimits, 'k--', 'LineWidth', 0.75);
    end
    fluxLegend{1 + k} = [num2str(targetFractions(k) * 100), '%: Age = ', num2str(round(targetAgesBP(k))), ' years BP'];
end
hold off
axis square
fig.Units = 'pixels';
fig.Position(3:4) = [figHeight, figHeight];
xlabel('Age [years BP]')
ylabel('Percentage of CH_4 flux younger than age')
legend(fluxLegend, 'Location', 'northwest')
exportgraphics(gca, fullfile(resultsDir, 'flux-age-distribution.png'));

%% Export flux-age distribution and modeled profiles to Excel
outFile = fullfile(resultsDir, 'flux-age-distribution.xlsx');
cols = {sedimentAgeClassesBP(:), 100 * cumulativeRelativeFlux(:), meshNodes, ...
        modeledProductionScaled(:, 1) ./ scaling, modeledMeanFm(:, 1), modeledConcentration(:, 1), ...
        methaneConcDepth, methaneConc, methaneFmDepth, methaneFm};
maxLen = max(cellfun(@length, cols));
for i = 1:numel(cols)
    cols{i}(end+1:maxLen) = NaN;
end
T = table(cols{1}, cols{2}, cols{3}, cols{4}, cols{5}, cols{6}, cols{7}, cols{8}, cols{9}, cols{10}, ...
    'VariableNames', {'year_BP', 'cumulative_percentage', 'z_m', ...
    'sim_prod_rate_mol_m3_s', 'sim_mean_Fm', 'sim_conc_mM', ...
    'z_conc_m', 'data_conc_mM', 'z_Fm_m', 'data_Fm'});
writetable(T, outFile);

%% Compute and display key model diagnostics
productionActual = productionMap' * productionCoeffs / scaling;
peakProduction = max(productionActual);
backgroundProduction = productionActual(end);

activeProduction = productionActual - backgroundProduction;
firstNearBackground = find(activeProduction < 1e-14, 1, 'first');
secondNearBackground = find(activeProduction < 1e-14, 2, 'first');
if numel(secondNearBackground) >= 2
    productionZoneThickness = meshNodes(secondNearBackground(2));
else
    productionZoneThickness = NaN;
end

[~, peakIdx] = max(productionActual);
peakProductionDepth = meshNodes(peakIdx);

interfaceConcentration = modeledConcentration(1);
interfaceGradient = (modeledConcentration(2) - modeledConcentration(1)) / (meshNodes(2) - meshNodes(1));
coeff0Phys = coeffFcn(0) / scaling;
totalFluxPhys = -coeff0Phys * interfaceGradient + qDarcy * interfaceConcentration;

dz = meshNodes(2) - meshNodes(1);
totalProductionPhys = sum(productionActual) * dz;

flux_mmol_m2_day = totalFluxPhys * 86400 * 1e3;
prod_mmol_m2_day = totalProductionPhys * 86400 * 1e3;
q_cm_yr = qDarcy * 100 * 86400 * 365;
peak_nmol_cm3_day = peakProduction * 86400 * 1e3;
bg_nmol_cm3_day = backgroundProduction * 86400 * 1e3;

fprintf('---------------------------------------\n');
fprintf('Model diagnostics\n');
fprintf('---------------------------------------\n');
fprintf('Upward methane flux at z=0 : %.3e mol m^-2 s^-1 (%.3f mmol m^-2 day^-1)\n', ...
    -totalFluxPhys, -flux_mmol_m2_day);
fprintf('Upward Darcy flux q        : %.3e m s^-1 (%.2f cm yr^-1)\n', ...
    -qDarcy, -q_cm_yr);
fprintf('Production zone thickness  : %.3f m\n', productionZoneThickness);
fprintf('Production peak depth      : %.3f m\n', peakProductionDepth);
fprintf('Production peak value      : %.3e mol m^-3 s^-1 (%.3f nmol cm^-3 day^-1)\n', ...
    peakProduction, peak_nmol_cm3_day);
fprintf('Background production      : %.3e mol m^-3 s^-1 (%.3f nmol cm^-3 day^-1)\n', ...
    backgroundProduction, bg_nmol_cm3_day);
fprintf('Total methane production   : %.3e mol m^-2 s^-1 (%.3f mmol m^-2 day^-1)\n', ...
    totalProductionPhys, prod_mmol_m2_day);
fprintf('---------------------------------------\n');

end
