function stimuli = facedist_importstimuli(subname)

% the root directory is 2 levels down (e.g. bidsdir/code/facedistid_analysis)
rootdir = fileparts(fileparts(fileparts(mfilename('fullpath'))));
if iscell(subname)
    % unpack AA5-style sub names (one cell per session)
    subname = subname{1}(1:6);
end
behaviourdir = fullfile(rootdir,'misc',subname);
assert(exist(behaviourdir,'dir')~=0,'could not find behaviourdir: %s',behaviourdir);

% process stimulusspace
stimpath = fullfile(behaviourdir,'stimuli_mri.mat');
ss = loadbetter(stimpath);
stimuli = struct('image',{ss.stimulus.imagealpha},'alpha',...
    {ss.stimulus.alpha},'shape',{ss.stimulus.shape},...
    'tex',{ss.stimulus.tex},'tl',{ss.stimulus.tl},...
    'azimuth',{ss.stimulus.azimuth});
% clean up alpha channel
newalpha = arrayfun(@(thisstim)im_isolatebiggest(thisstim.image>254,...
    'alpha',thisstim.alpha),stimuli,'uniformoutput',false);
[stimuli.alpha] = deal(newalpha{:});

% NEW - also crop the image - first by masking images 1:12 and 13:24
% separately (to a common size), then applying said mask to the images.
%
% apply a coloured background around each according to ecc colours
% (here we may need to consider resortind to get it right)
leftmask = matmean(stimuli(1:12).alpha);
rightmask = matmean(stimuli(13:24).alpha);
leftmask(:,end-1:end) = 0;
% I don't think there is a similar problem on this side but just to be
% safe...
rightmask(:,1:2) = 0;
imsize = size(rightmask);

% now trim the columns
leftcol = find(sum(leftmask,1)>0);
leftrange = range([leftcol(1) leftcol(end)]);
rightcol = find(sum(rightmask,1)>0);
rightrange = range([rightcol(1) rightcol(end)]);
ntouse = ceil(max([leftrange rightrange])*1.1);

% take the first ntouse for left, and the last ntouse for right.
leftind = 1:ntouse;
rightind = leftind + (imsize(2)-ntouse);

scalef = [1.05 1.2];
finalsize = ceil([imsize(1) ntouse] .* scalef);
oval = imresize(double(makecircle(1000)),finalsize);
% smooth slightly to make the border look more continuous
oval = imfilter(oval,fspecial('gaussian',20,2),'same');
oval = repmat(oval,[1 1 3]);

% 3D settings
proptrim = .5;
nfaces = 500;

% go through and apply hard alpha blend of face image with solid
% colour background
bg = ones([imsize(1) ntouse 3]);
inds = [repmat({leftind},[12 1]); repmat({rightind},[12 1])];
colors = repmat(facedist_colors('ecc'),[8 1]);
for n = 1:24
    ind = inds{n};
    % trim the edge and set type
    thisim = im2double(stimuli(n).image(:,ind,:));
    thisalpha = stimuli(n).alpha(:,ind,:);
    % hang on to the raw image and alpha values for image processing
    % analyses
    stimuli(n).rawimage = thisim;
    stimuli(n).rawalpha = thisalpha;
    % make rgb
    thisalpha = repmat(thisalpha,[1 1 3]);
    thisim = repmat(thisim,[1 1 3]);
    % hard blend in colored background
    thisbg = bsxfun(@times,bg,reshape(colors(n,:),[1 1 3]));
    thisim = (thisim .* thisalpha) + (thisbg .* (1-thisalpha));
    % upsize to final (basically add some blank space for the oval to show
    thisim = im_canvas(thisim,finalsize,colors(n,:));
    thisalpha = im_canvas(thisalpha,finalsize,0);
    % finally, hard blend in a white background in the oval
    thisim = (thisim .* oval) + (1-oval);
    % and now the alpha channel is simply the oval.
    stimuli(n).image = thisim;
    stimuli(n).alpha = oval(:,:,1);
    % make sure 0 is 0
    stimuli(n).image(stimuli(n).image<0) = 0;
    stimuli(n).alpha(stimuli(n).alpha<0) = 0;
    stimuli(n).color = colors(n,:);

    % also handle texture details
    % first trim posterior
    ymin = min(stimuli(n).shape(:,2));
    cutoff = range(stimuli(n).shape(:,2))*proptrim+ymin;
    % now the question is, can we just drop bad indices?
    bad = stimuli(n).shape(:,2)<cutoff;
    % the key challenge is that tl includes indices into shape and tex
    % rows. So when we delete rows we need to adjust the indexing to
    % match.
    % find the bad cases - cases where y is < cutoff
    badtl = stimuli(n).shape(stimuli(n).tl(:,2),2) < cutoff;
    % so the rowindices used to be 1:n, not we want to count 1:n but
    % skipping bad ones
    ind = ones(size(stimuli(n).shape,1),1);
    ind(bad) = 0;
    ind = cumsum(ind);
    % so now the trick is to update all the indices in tl to match.
    % basically, any index 1:n in the original should be replaced with the
    % vector here. 
    % so
    stimuli(n).tl(:,1) = ind(stimuli(n).tl(:,1));
    stimuli(n).tl(:,2) = ind(stimuli(n).tl(:,2));
    stimuli(n).tl(:,3) = ind(stimuli(n).tl(:,3));
    % and now removing the bad rows should be safe
    stimuli(n).tl(badtl,:) = [];
    % find points in stimuli(n).tl that include y values < cutoff
    stimuli(n).shape(bad,:) = [];
    stimuli(n).tex(bad,:) = [];
    
    % center on 0 on all axes
    stimuli(n).shape = bsxfun(@minus,stimuli(n).shape,...
        mean(stimuli(n).shape));
    % scale to 0-1 range
    maxs = max(stimuli(n).shape(:));
    stimuli(n).shape = stimuli(n).shape / maxs;

    % to rotate the mesh we would need to apply some kind of transformation
    % matrix to the coordinates. Let's revisit

    % DEBUG
    %F = figure;
    %h = trimesh(stimuli(n).tl,stimuli(n).shape(:,1),...
        %stimuli(n).shape(:,2),stimuli(n).shape(:,3),...
        %'EdgeColor','none','FaceVertexCData',stimuli(n).tex,...
        %'facecolor','interp','facelighting','phong');
    % this seems to work so far.
    %hs = addellipsoid([1.15 2 1.2],[0 -.73 .2],...
        %range(stimuli(n).shape,1),stimuli(n).color,nfaces)
    %set(gca,'projection','perspective','view',[150 0],'dataaspectratio',...
        %[1 1 1]);
    %keyboard;
end
