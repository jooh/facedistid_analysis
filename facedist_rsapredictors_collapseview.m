function noviewrdms = facedist_rsapredictors_collapseview(rdms)

% now get rid of the view-specific ones
hits = strfindcell({rdms.name},'view');
rdms(hits) = [];

% and for the rest, average over the two view matrices (this doesn't really
% change anything if they're the same of course)

for r = 1:length(rdms)
    rdmstack = {rdms(r).RDM(1:12,1:12),rdms(r).RDM(13:24,13:24)};
    noviewrdms(r) = rdms(r);
    noviewrdms(r).RDM = matmean(rdmstack{:});
end
