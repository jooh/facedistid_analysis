% generate a sample with a 2D chi distribution. df (which can be
% continuous!) interpolates from gaussian (1) to more donut-like values.
%
% xy = facedist_chi2sample(nunit,df)
function xy = facedist_chi2sample(nunit,df)

% so according to niko, you can just sample from the inv cdf
radii = chi2rnd(df,nunit,1);
[x,y] = pol2cart(rand([nunit,1])*2*pi,sqrt(radii));
xy = [x y];
