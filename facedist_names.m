function names = facedist_names(varargin)

names = {};

for n = 1:nargin
    switch varargin{n}
        case 'delta_eccentricity'
            names(end+1) = {'\Deltaeccentricity'};
        case 'delta_direction'
            names(end+1) = {'\Deltadirection'};
        case 'full'
            names(end+1) = {'face space distance'};
        case {'facepairs','fit: facepairs'}
            names(end+1) = {'perceptual judgments'};
        case {'GWP','gwp','gwpgrid'}
            names(end+1) = {'gabor filter'};
        case {'gwpgrid_noav'}
            names(end+1) = {'gabor filter (a=0)'};
        case 'rep_lower'
            names(end+1) = {'noise ceiling (lower bound)'};
        case 'rep_upper'
            names(end+1) = {'noise ceiling (upper bound)'};
        case 'vid_corrdist'
            names(end+1) = {'pixelwise correlation'};
        case {'FFA','fit: FFA'}
            names(end+1) = {'fusiform face area'};
        case {'EVC','fit: EVC'}
            names(end+1) = {'early visual cortex'};
        case {'OFA','fit: OFA'}
            names(end+1) = {'occipital face area'};
        case {'PPA','fit: PPA'}
            names(end+1) = {'parahippocampal place area'};
        case {'TOS','fit: TOS'}
            names(end+1) = {'transverse occipital sulcus'};
        case {'RSC','fit: RSC'}
            names(end+1) = {'retrosplenial cortex'};
        case 'ex_gaussian'
            names(end+1) = {'gaussian exemplar'};
        case 'ex_gaussian_noav'
            names(end+1) = {'gaussian exemplar (a=0)'};
        case 'ex_negativegaussian'
            names(end+1) = {'negative gaussian exemplar'};
        case 'ex_negativegaussian_noav'
            names(end+1) = {'negative gaussian exemplar (a=0)'};
        case 'ramp'
            names(end+1) = {'sigmoidal ramp tuning'};
        case 'ramp_noav'
            names(end+1) = {'sigmoidal ramp tuning (a=0)'};
        case 'view_across_sub-sub_m6.0'
            names(end+1) = {'view across: sub-caricature'};
        case 'view_within_sub-sub_m6.0'
            names(end+1) = {'view within: sub-caricature'};
        case 'sub-sub_m6.0'
            names(end+1) = {'sub-caricature'};
        case 'view_across_typical-typical_m19.9'
            names(end+1) = {'view across: typical'};
        case 'view_within_typical-typical_m19.9'
            names(end+1) = {'view within: typical'};
        case 'typical-typical_m19.9'
            names(end+1) = {'typical'};
        case 'view_across_caricature-caricature_m33.8'
            names(end+1) = {'view across: caricature'};
        case 'view_within_caricature-caricature_m33.8'
            names(end+1) = {'view within: caricature'};
        case 'caricature-caricature_m33.8'
            names(end+1) = {'caricature'};
        case 'view_across_slope_absecc'
            names(end+1) = {'view across: eccentricity slope'};
        case 'view_within_slope_absecc'
            names(end+1) = {'view within: eccentricity slope'};
        case 'slope_absecc'
            names(end+1) = {'eccentricity slope'};
        otherwise
            % just return the same name
            names(end+1) = varargin(n);
    end
end
