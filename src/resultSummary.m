function t = resultSummary(r)
%RESULTSUMMARY One-row table; all fluxes are per unit sediment surface area.
d=r.diagnostics; s=struct();
s.q_m_s=r.qMagnitude; s.q_cm_year=r.qMagnitude*100*86400*365;
s.k_loss_day_inv=r.kLoss_per_day;
s.z_peak_m=r.productionFit.shape(1); s.w_tail_m=r.productionFit.shape(2);
s.P0_mol_m3_bulk_s=r.productionFit.amplitudesScaled(1)/r.config.scaling;
s.Pbg_mol_m3_bulk_s=r.productionFit.amplitudesScaled(2)/r.config.scaling;
s.k_sw_m_s=r.kSw; s.z_loss_m=r.config.lossDepth;
s.E_C=r.E_C; s.E_F=r.E_F;
s.RMSE_C_mM=d.rmseConcentration; s.RMSE_Fm=d.rmseFm;
s.J_gross_mmol_m2_day=d.grossProduction*86400*1000;
s.J_loss_mmol_m2_day=d.loss*86400*1000;
s.J_export_mmol_m2_day=d.upwardFlux*86400*1000;
s.J_export_diffusive_mmol_m2_day=d.diffusiveUpwardFlux*86400*1000;
s.J_export_advective_mmol_m2_day=d.advectiveUpwardFlux*86400*1000;
s.loss_fraction_of_gross=d.lossFraction;
s.legacy_fraction_of_export=d.legacyFraction;
s.legacy_export_mmol_m2_day=d.legacyExport*86400*1000;
s.median_source_apparent_age_BP=d.ageQuantiles(3);
s.surface_Fm=d.surfaceFm;
t=struct2table(s);
end
