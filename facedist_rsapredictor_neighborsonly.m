% Return RDMs for plotting a dissimilarity grid with plotmdslines
%
% [dirrdm eccrdm] = facedist_rsapredictor_neighborsonly(ncon)
function [dirrdm eccrdm] = facedist_rsapredictor_neighborsonly(ncon)
% make an RDM where only the neighbouring directions are non-NaN (input for
% plotmdslines)
eccrdm = zeros([12 12]);
for x = 1:3
    eccrdm(x,x+1) = 1;
    eccrdm(x+1,x) = 1;
end
for y = 2:3
    eccrdm = eccrdm + y*circshift(eccrdm,[4 4]);
end
eccrdm(eccrdm==0) = NaN;
eccrdm = zerodiagonal(eccrdm);

% similarly, we want an RDM where only the neighbouring eccentricities are
% non-NaN
baseoff = zeros(4);
basemat = eye(4);
dirrdm = [baseoff basemat baseoff; basemat baseoff basemat; ...
    baseoff basemat baseoff];
dirrdm(dirrdm==0) = NaN;
dirrdm = zerodiagonal(dirrdm);

if ncon == 24
    % upcast to two views
    fullmat = zerodiagonal(NaN([24 24 2]));
    fulldir = fullmat;
    fulldir(1:12,1:12,1) = dirrdm;
    fulldir(13:24,13:24,2) = dirrdm;
    dirrdm = fulldir;
    fullecc = fullmat;
    fullecc(1:12,1:12,1) = eccrdm;
    fullecc(13:24,13:24,2) = eccrdm;
    eccrdm = fullecc;
end
