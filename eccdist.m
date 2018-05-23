% an eccentricity (from origin) distance function to be used with
% pdist(x,@eccdist).
% D = eccdist(XI,XJ);
function D = eccdist(XI,XJ);

XIdist = distorigin(XI);
XJdist = distorigin(XJ);

D = abs(XJdist-XIdist);
