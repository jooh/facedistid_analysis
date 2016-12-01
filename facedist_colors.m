% standardised colour schemes for facedistid.
%
% c = facedist_colors(varargin)
function c = facedist_colors(varargin)

pastels = cmap_pastel - .2;

% color gradient for HMAX. Hm.
hmax = colorScale([0 .1 .5; .2 .5 1],4);

c = [];
eccs = circshift(cmap_redyellow(3),[0 2]) * .95;
% something toward green
v = [.2 .9 .7];
for n = 1:nargin
    switch varargin{n}
        case {'EVC','fit: EVC'}
            c = [c; pastels(11,:)];
        case {'FFA','fit: FFA'}
            % c = [c; pastels(4,:)];
            c = [c; pastels(3,:)];
        case {'OFA','fit: OFA'}
            c = [c; pastels(5,:)];
        case {'PPA','fit: PPA'}
            c = [c; pastels(6,:)];
        case {'TOS','fit: TOS'}
            c = [c; pastels(10,:)];
        case {'facepairs','fit: facepairs'}
            c = [c; pastels(4,:)];
        case 'fmri'
            c = [c; pastels(2,:)];
        case 'beh'
            c = [c; pastels(2,:)];
        case {'ramp','ex_gaussian','full','ex_negativegaussian','ramp_noav','ex_gaussian_noav','ex_negativegaussian_noav'}
            c = [c; pastels(1,:)];
        case 'full'
            c = [c; [.4 .4 .4]];
        case 'hmax_c1'
            c = [c; hmax(1,:)];
        case 'hmax_c2'
            c = [c; hmax(2,:)];
        case 'hmax_c2b'
            c = [c; hmax(3,:)];
        case 'hmax_c3'
            c = [c; hmax(4,:)];
        case {'gwpgrid','gwpgrid_noav','vid_corrdist'}
            c = [c; pastels(5,:)];
        case 'rsaglm'
            c = [c; [.9 0 0.1]];
        case {'sub','sub-sub'}
            c = [c; eccs(1,:)];
        case {'typical','typical-typical'}
            c = [c; eccs(2,:)];
        case {'caricature','caricature-caricature'}
            c = [c; eccs(3,:)];
        case 'sub-typical'
            c = [c; mean(eccs(1:2,:))];
        case 'typical-caricature'
            c = [c; mean(eccs(2:3,:))];
        case 'sub-caricature'
            c = [c; mean(eccs([1 3],:))];
        case 'ecc'
            c = [c; eccs];
        case 'shadecolor'
            c = [c; [.8 .8 .8]];
        case 'fitlinecolor'
            c = [c; [0 0 0]];
        case {'eccentricity','eccentricity*direction'}
            c = [c; mean(eccs)];
        case 'direction'
            c = [c; 1 1 1];
        case {'view','view*direction'}
            c = [c; v];
        case {'view*eccentricity','view*eccentricity*direction'}
            c = [c; mean([v; mean(eccs)])];
        case 'reproducibility'
            c = [c; [.8 .8 .8]];
        case 'singles'
            c = [c; [.5 .5 .5]];
        otherwise
            error('no color for: %s',varargin{n});
    end
end
