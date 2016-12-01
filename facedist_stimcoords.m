% obtain the unique radii and polar coordinates for facedistid
%
% [phi,rad] = facedist_stimcoords(assumesymmetry)
function [phi,rad] = facedist_stimcoords(assumesymmetry)

if assumesymmetry
    phi = linspace(0,2*pi,7);
else
    phi = linspace(0,pi,4);
end

% rotate the whole thing
phi = rem(phi+1.5*pi,2*pi);

rad = [4.232, 14.107, 23.981];


