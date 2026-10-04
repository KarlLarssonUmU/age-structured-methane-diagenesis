function data = loadModelData(cfg)
%LOADMODELDATA Read the explicit, unit-labelled main-model CSV inputs.
validateConfig(cfg);
names={'sediment_porosity.csv','methane_concentration.csv', ...
    'sediment_source_reservoir_corrected.csv','methane_fm_observations.csv'};
data.inputFiles=cellfun(@(n) fullfile(cfg.repoDir,'data',n),names,'UniformOutput',false);
p=readtable(data.inputFiles{1}); c=readtable(data.inputFiles{2});
s=readtable(data.inputFiles{3}); f=readtable(data.inputFiles{4});
data.porosityDepth=p.depth_m; data.porosity=p.porosity_fraction;
data.concentrationDepth=c.depth_m; data.concentration=c.CH4_mM;
data.sedimentDepth=s.depth_m; data.sedimentAge=s.age_BP;
data.sedimentFm=s.source_Fm;
data.observedFmDepth=f.depth_cm/100; data.observedFm=f.Fm;
data.observedFmSd=f.Fm_sd;
data.fmDepth=data.observedFmDepth; data.fm=data.observedFm;
for name={'porosityDepth','porosity','concentrationDepth','concentration', ...
        'sedimentDepth','sedimentAge','sedimentFm','fmDepth','fm','observedFmSd'}
    v=data.(name{1});
    assert(iscolumn(v) && ~isempty(v) && all(isfinite(v)), ...
        'methane:InputData','Nonfinite or empty input column: %s',name{1});
end
assert(all(diff(data.porosityDepth)>0) && all(diff(data.sedimentDepth)>0) && ...
    all(diff(data.concentrationDepth)>0) && all(diff(data.fmDepth)>=0), ...
    'methane:InputData','Depths must be ordered; profile depths must be unique.');
assert(data.sedimentDepth(1)<=0 && data.sedimentDepth(end)>=cfg.bottomDepth && ...
    data.porosityDepth(1)<=0 && data.porosityDepth(end)>=cfg.bottomDepth, ...
    'methane:InputData','Porosity and source curves must cover the whole domain.');
assert(all(data.porosity>0 & data.porosity<=1) && ...
    all(data.concentration>0) && all(data.fm>0) && all(data.observedFmSd>0));
assert(all(diff(data.sedimentAge)>=0) && all(data.sedimentFm>0));
assert(max(abs(data.sedimentFm-exp(-data.sedimentAge/8033)))<1e-10, ...
    'methane:SourceFm','Source Fm must match the supplied adjusted apparent ages.');
assert(all(data.concentrationDepth>=0 & data.concentrationDepth<=cfg.bottomDepth) && ...
    all(data.fmDepth>=0 & data.fmDepth<=cfg.bottomDepth));
assert(cfg.concentrationFitPoints<=numel(data.concentration));
data.surfaceLine=polyfit(data.concentrationDepth(1:cfg.concentrationFitPoints), ...
    data.concentration(1:cfg.concentrationFitPoints),1);
assert(all(data.surfaceLine>0),'methane:SurfaceFit','Invalid near-surface regression.');
end
