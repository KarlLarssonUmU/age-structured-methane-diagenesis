function B = productionBasis(z,family,shape)
%PRODUCTIONBASIS Columns multiply [P0; Pbg], the two fitted amplitudes.
assert(strcmp(family,'bump'),'methane:ProductionFamily', ...
    'This release implements the constrained smooth bump.');
B=productionShape(z,shape(1),shape(2))';
end
