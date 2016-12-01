function crossrdms = facedist_rsapredictors_crossview(rdms)

% now get rid of the view-specific ones
hits = strfindcell({rdms.name},'view');
rdms(hits) = [];

% and pick the off quadrant (setting the diagonal to 0 as we go)
for r = 1:length(rdms)
    % average the two quadrants (to get A>B and B>A mean)
    rdmstack = {rdms(r).RDM(13:24,1:12),rdms(r).RDM(1:12,13:24)};
    crossrdms(r) = rdms(r);
    crossrdms(r).RDM = matmean(rdmstack{:});
    % make sure diagonal is 0
    crossrdms(r).RDM(logical(eye(12))) = 0;
end
