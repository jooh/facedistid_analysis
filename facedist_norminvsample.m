function xy = facedist_norminvsample(nunit,v)

% we have a uniform distribution of directions
thetas = rand(nunit,1) * 2 * pi;

% now we want to generate eccentricities from uniform and chi2
% distributions. The insight from Niko is that you can use the inverse
% distributions for this.
randecc = rand(nunit,1);
% get chi2 eccentricities for these values
if v>0
    chiecc = chi2inv(randecc,2);
else
    % need to flip here to match the flip of the distribution later.
    % Otherwise the intermediate morphs tend to have narrow distance
    % distributions (since you are morphing between anti-correlated values)
    chiecc = chi2inv(1-randecc,2);
end
% then uniform eccentricities in a similar range
% this value determines how much the distributions get truncated at the
% ends, and how far out the inversion happens for the gaussian. So there is
% a trade-off between not truncating the distribution too badly, and how
% compressed the inverted gaussian ends up.
criticalchi = chi2inv(.99,2);
uniecc = unifinv(randecc,0,criticalchi);

% handle negative gaussian case
if v < 0
    chiecc = criticalchi - chiecc;
end

% interpolate the two vectors
finalecc = morph(chiecc,uniecc,abs(v));
% rescore negatives
finalecc(finalecc<0) = 0;
% truncate positive values
finalecc(finalecc>criticalchi) = criticalchi;

% scale so that the final distances will end at the caricatures
finalecc = sqrt(finalecc) *  (2.4042/sqrt(criticalchi));

% generate the cartesian coordinates
[x,y] = pol2cart(thetas,finalecc);
xy = [x(:),y(:)];

