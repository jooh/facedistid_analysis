% reference coordinates for facedistid MDS plots. With 1 view (ncon=12) we
% plot the spokes extending out to the right, with 2 views (ncon=24) we
% plot the second with with spokes extending out to the left. When used
% with procrustes this ensures that all plots have roughly the same
% orientation.
%
% xy = facedist_getrefxy(ncon)
function xy = facedist_getrefxy(ncon)

% coordinates with the face space norm roughly at [0 0] and spokes
% extending towards the right
xy = [...
    0.0006    4.2320;
    3.6655    2.1160;
    3.6655   -2.1160;
    0.0006   -4.2320;
    0.0010   14.1067;
    12.2173    7.0532;
    12.2173   -7.0532;
    0.0010  -14.1067;
    0.0014   23.9814;
    20.7691   11.9904;
    20.7691  -11.9904;
    0.0014  -23.9814];
switch ncon
    case 12
        % done
    case 24
        xy = [xy; xy];
    otherwise
        error('unknown ncon: %d',ncon);
end
