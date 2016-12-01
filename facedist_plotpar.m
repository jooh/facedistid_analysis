function par = facedist_plotpar

par.pthresh_singles = .05;
par.abseccfsize = [10 6] * 1.3;
par.dirfsize = par.abseccfsize .* [1 1];
par.eccdirfsize = par.abseccfsize * sqrt(1.6) .* [.7 1];
par.multifsize = par.abseccfsize * 1.2;
par.rsaglmfsize = [5 8];
par.fitmdsfsize = [8 8] * 1.7;
par.fitratiofsize = [10 7.5];
par.fitunisize = [6 6];
par.axscale = 1/3;
par.basecolors = facedist_colors('ecc');
par.shadecolor = [.7 .7 .7];
par.notsigcolor = [.6 .6 .6];
par.fitlinecolor = [0 0 0];
par.fitlinewidth = 1;
par.chancelinewidth = .5;
par.chancelinestyle = ':';
par.chancelinecolor = [0 0 0];
par.barwidth = .5;
par.errorbarlinewidth = .5;
par.markersize = 4;
par.markertype = 'o';
par.gridlines =[4 8 12 16 20];
par.ticklength = [.02 .02];
par.subxscale = par.barwidth .* [-.25 .25];
par.eccxoff = par.barwidth .* [-.4 0 .4];
par.viewxoff = [-.2 +.2];
par.absecclabels = {'sub','typical','caricature'};
par.abseccxlabel = 'eccentricity';
par.rsaglmxlabel = 'view';
par.dirxlabel = '\Deltadirection \circ';
par.dirxlabel_glm = 'direction \circ';
par.dirlabels = {'0','60','120','180'};
par.ylabel = {'face exemplar information',...
    '(cross-validated mahalanobis distance)'};
par.ylabel_rsaglm = {'face space contribution',...
    '(parameter estimate)'};
par.ylabel_fitratio = {'face space warp index'};
par.ylabel_facepairs = {'rated dissimilarity',...
    '(% judged more different)'};
par.ylabel_glm = {'univariate fMRI response (a.u.)',...
    '(parameter estimate)'};
par.outlinecolor = 'm';
par.dochanceshade = false;
