classdef FaceFigure < hgsetget & dynamicprops
    % Display a face mesh based on shape/texture information from a
    % BaselSpace model Returns an object with handles and methods for
    % on-the-fly changes to view/lighting and export to images. Typically
    % invoked through the BaselSpace.showface method. varargin allows
    % setting of ANY of the class properties - caveat emptor.
    %
    % f = FaceFigure(tl,shape,tex,segbin,[varargin])
    %
    % FaceFigure Properties:
    % shapecoef = []; % original PC loading used to construct face shape
    % texcoef = []; % original PC locading used to construct face tex
    % tl = []; % tl field from BaselSpace
    % shape = []; % xyz vertices
    % tex = []; % rgb colour for each vertex
    % fighand=[]; % figure handle
    % axhand=[]; % axis handle
    % mhand=[]; % mesh handle
    % azimuth=0; % face horizontal view
    % elevation=0; % face vertical view
    % figsize = [640 640]; % figure size in pixels
    % bgcolor = ones(1,3)*.5; % background colour
    % lights = struct('hand',{},'azimuth',{},'elevation',{},...
        % 'style',{},'color',{},'centeroncam',{});
    % ambientstrength = .5; % ambient light intensity
    % diffusestrength = .5; % affects light source intensity
    % specularstrength = .1; % specular intensity from light sources
    % specularexponent = 1;
    % zoom; % camva setting (customise with scalar in deg)
    % centercamoneyes = true; % optional setting to center camera on eyes
    % trimfaceprop = 0; % optionally call removeface on init
    % exportwithannulus = 1; % include an annulus in exported ims
    % visible = 'on'; % fig visible on initialise
    % image % exported image form of figure (new on every updateface)
    % alpha % alpha annulus (new on every openfigure)
    % imagealpha % image with 'hard' alpha blend
    % scramble % scrambled version of image (new on every openfigure)
    % scrambletarget = 'image' % target for scrambling
    % grayscale = 0; % flag to detect grayscale mode (no colour
    % channels)
    %
    % FaceFigure Methods:
    % turnface - rotate the face view
    % addlight - add a light source
    % turnlight - change an existing light source's position
    % clearlight - remove some or all current lights
    % defaultambientlight - standard ambient light settings
    % hardambientlight - 'spotlight' effect with lots of shadows
    % onlyambientlight - kill non-ambient light sources
    % changetexture - switch texture between original and segmentation
    % changezoom - zoom in or out in degrees
    % removeposterior - remove back of the mesh from view
    % exportasmatrix - export current figure to image matrix
    % exportasimage - write current figure to a PNG file

    properties
        tl = []; % tl field from BaselSpace
        shape = []; % shape mesh
        tex = []; % texture mesh
        shapecoef = []; % shape coefficient
        texcoef = []; % tex coefficient
        fighand=[]; % figure handle
        axhand=[]; % axis handle
        mhand=[]; % mesh handle
        azimuth=0; % face horizontal view
        elevation=0; % face vertical view
        figsize = [640 640]; % figure size in pixels
        bgcolor = ones(1,3)*.5; % background colour
        lights = struct('hand',{},'azimuth',{},'elevation',{},...
            'style',{},'color',{},'centeroncam',{});
        ambientstrength = .5; % ambient light intensity
        diffusestrength = .5; % affects light source intensity
        specularstrength = .1; % specular intensity from light sources
        specularexponent = 1;
        zoom; % camva setting (customise with scalar in deg)
        centercamoneyes = true; % optional setting to center camera on eyes
        trimfaceprop = 0; % optionally call removeface on init
        exportwithannulus = 1; % include an annulus in exported ims
        visible = 'on'; % fig visible on initialise
        image = []; % RGB image of current figure
        alpha = []; % alpha annulus layer for image
        imagealpha = []; % hard alpha blend of image
        scramble = []; % scrambled image
        scrambletarget = 'image'; % scramble target
        alphamask = true; % use alpha mask?
        imageautotrim = true; % trim the image output
        imagelims = []; % store limits for image trimming
        texorg; % store original face texture 
        segbin; % keep track of where different image parts are
        videoframes = []; % store video (output from rotateface)
        grayscale = 0; % flag to detect grayscale mode 
        eyesxyz; % center of mass of eye region for centercamoneyes
        defxyz; % default camera center for centercamoneyes==0
    end

    properties(Access = private)
        azioffset=180; % azimuth offset to center view on 0 deg
        lightsdefault = struct('hand',NaN,'azimuth',0,'elevation',0,...
            'style','local','color',[1 1 1],...
            'centeroncam',false);
        seglab2ind = struct('nose',1,'eyes',2,'mouth',3,'external',4);
        yoffset = .09;
        yoffsetalt = -.05;
        xoffset = 0; %-.004;
        fh = .5;
        fw = .35;
        cuttop = .15;
        zoomfactor = 10;
    end

    methods
        function f = FaceFigure(tl,shape,tex,segbin,varargin)
        % Create a face figure. For a complete list of possible extra
        % arguments, see get(f) after creating instance or help FaceFigure.
        % f = FaceFigure(shape,tex,varargin)
            if nargin==0
                % support default initialisation mode for object arrays
                return
            end
            f = varargs2structfields(varargin,f); 
            % Parse string inputs for visibility
            if ~ischar(f.visible)
                if f.visible
                    f.visible = 'on';
                else
                    f.visible = 'off';
                end
            end
            f.tl = tl;
            f.shape = shape;
            f.tex = tex;
            f.segbin = segbin;
            if isempty(f.fighand)
                % find a unique handle
                f.fighand = uniquehandle;
            elseif ishandle(f.fighand)
                % make sure we don't start drawing into an existing figure
                % since this rarely works well
                close(f.fighand);
            end
            f.openfigure;
        end

        function changezoom(self,newzoom)
        % Change the zoom (ie scale) in degrees.
            self.zoom = newzoom;
            self.updateface;
        end

        function addlight(self,azimuth,elevation,varargin)
        % Add a light source. Azimuth and elevation are relative to the
        % current head. For possible varargins and defaults, see
        % fieldnames(self.lights)
        % addlight(self,azimuth,elevation,varargin)
            lind = length(self.lights)+1;
            %fprintf('(FaceFigure) adding light %d to figure %d\n',...
                %lind,self.fighand)
            % Initialise with default light settings
            % (extra check to support the case where lights==[]
            if lind==1
                self.lights = self.lightsdefault;
            else
                self.lights(lind) = self.lightsdefault;
            end
            % And update with inputs
            self.lights(lind).azimuth = azimuth;
            self.lights(lind).elevation = elevation;
            self.lights(lind) = varargs2structfields(varargin,self.lights(lind));
            self.updateface;
        end

        function turnlight(self,azimuth,elevation,lind)
        % Change azimuth and elevation of a currently displayed light in
        % degrees relative to current view. lind is an index into
        % self.lights (default 1)
        % turnlight(self,azimuth,elevation,[lind])
            if ieNotDefined('lind')
                lind = 1;
            end
            self.lights(lind).azimuth = azimuth;
            self.lights(lind).elevation = elevation;
            self.updateface;
        end

        function turnface(self,azimuth,elevation)
        % Change the azimuth and elevation of the currently displayed face
        % in degrees. Unlike the built-in Matlab view function, here
        % angles=[0 0] produces a frontal view. 
        % turnface(self,azimuth,elevation)
            self.azimuth = azimuth;
            self.elevation = elevation;
            self.updateface;
        end

        function frames = rotateface(self,azilims,elelims,nframes,target)
        % Change the azimuth/elevation of a face in nframes linear
        % increments, return the result as a 4D r g b frame matrix. The
        % result is also stored under ff.videoframes
        % azilims: undefined, scalar (ie don't change) or [1 2] vector
        % defining star and end azimuth.
        % elelims: undefined, scalar (ie don't change) or [1 2] vector
        % defining star and end elevation.
        % nframes: number of steps (including start and end).
        % target: default image, but imagealpha may also be interesting.
        % frames = rotateface(self,[azilims],[elelims],nframes,[target])
            if ieNotDefined('azilims')
                % no azimuth change - centered
                azilims = [0 0];
            elseif length(azilims)==1
                % no azimuth change - fixed azimuth
                azilims = repmat(azilims,[1 2]);
            elseif length(azilims)>2
                error('azilims must be undefined, scalar, or length 2')
            end
            azisteps = linspace(azilims(1),azilims(2),nframes);
            if ieNotDefined('elelims')
                elelims = [0 0];
            elseif length(elelims)==1
                elelims = repmat(elelims,[1 2]);
            elseif length(elelims)>2
                error('elelims must be undefined, scalar, or length 2')
            end
            if ieNotDefined('target')
                target = 'image';
            end
            elesteps = linspace(elelims(1),elelims(2),nframes);
            % insert singleton 3rd dim for BW targets
            switch ndims(self.(target))
                case 2
                    vidsize = [size(self.(target)) 1];
                case 3
                    vidsize = size(self.(target));
                otherwise
                    error('unrecognised target dimensionality')
            end
            self.videoframes = zeros([vidsize nframes],'uint8');
            for n = 1:nframes
                self.turnface(azisteps(n),elesteps(n));
                self.videoframes(:,:,:,n) = self.(target);
            end
            if nargout > 0
                frames = self.videoframes;
            end
        end

        function clearlight(self,lind)
        % Remove light source from plot. lind gives scalar indices into
        % self.lights (default all)
        % clearlight(self,[lind])
            if isempty(self.lights)
                % Nothing to clear
                return
            end
            if ieNotDefined('lind')
                lind = 1:length(self.lights);
            end
            for l = 1:length(lind)
                delete(self.lights(lind(l)).hand);
            end
            % update record
            self.lights(lind) = [];
        end

        function defaultambientlight(self)
        % Set ambient light levels of current figure to the standard. Bit
        % of shadow but not too harsh.
        % defaultambientlight()
            self.ambientstrength = .5;
            self.diffusestrength = .5;
            self.specularstrength = .08;
            self.updateface;
        end

        function hardambientlight(self)
        % Set ambient light levels of current figure to cast hard shadows.
        % hardambientlight()
            self.ambientstrength = .3;
            self.diffusestrength = .7;
            self.specularstrength = .08;
            self.updateface;
        end

        function onlyambientlight(self)
        % Kill all light sources in current figure, set ambient to one. No
        % shadows at all.  
        % onlyambientlight()
            self.clearlight;
            self.ambientstrength = 1;
            self.diffusestrength = 0;
            self.specularstrength = 0;
        end

        function changetexture(self,textype)
        % Switch the vertex colormap to show the segmentation or the
        % original colour map. textype is 'segmentation', 'original',
        % 'nose', 'eyes', 'mouth', or 'external'.
        % changetexture(self,textype)
            [nvert,nseg] = size(self.segbin);
            segrgb = ones(nvert,3)*.5;
            switch lower(textype)
                case 'segmentation'
                    % Resample segment matrix to 3 columns (RGB)
                    segs = resample(double(self.segbin)',3,nseg)';
                    segrgb(segs>0) = segs(segs>0);
                case 'original'
                    segrgb = self.texorg;
                case 'grayscale'
                    segrgb = repmat(mean(self.tex,2),[1 3]);
                otherwise
                    segi = self.seglab2ind.(lower(textype));
                    segs = repmat(self.segbin(:,segi),[1 3]);
                    segrgb(segs>0) = segs(segs>0);
            end
            self.tex = segrgb;
            self.updateface;
        end

        function removeposterior(self,prop)
        % Remove the back of the mesh from view (effectively trimming ears,
        % neck and cheek). prop defines the cutoff point as a proportion of
        % the depth of the face (0 removes nothing, 1 everything).
        % removeposterior(self,prop)
            verts = get(self.mhand,'vertices');
            ymin = min(verts(:,2));
            cutoff = range(verts(:,2))*prop+ymin;
            set(self.mhand,'facealpha','flat','facevertexalphadata',...
                double(verts(:,2)>cutoff));
        end

        function [im ann] = exportasmatrix(self)
        % Export the current figure to a [x y 3] uint8 matrix. Faster than
        % exporting to a file with print and im_read'ing the result, but
        % resolution is limited by figsize. See also exportasimage.
        % [im alpha] = exportasmatrix()
            orgppm = get(self.fighand,'paperpositionmode');
            set(self.fighand,'paperpositionmode','auto');
            im = hardcopy(self.fighand,'-Dopengl','-r0');
            if self.exportwithannulus && nargout==2
                ann = self.addannulus(im);
            elseif nargout==2
                error('asked for alpha output but exportwithannulus is off')
            end
            set(self.fighand,'paperpositionmode',orgppm);
            if self.grayscale
                im = uint8(mean(im,3));
            end
        end

        function exportasimage(self,filename,resolution)
        % Export the current figure to a png. Resolution defaults to 300
        % (dpi). See also exportasmatrix.
        % exportasimage(self,filename,[resolution])
            if ieNotDefined('resolution')
                resolution = 300;
            end
            orgppm = get(self.fighand,'paperpositionmode');
            set(self.fighand,'paperpositionmode','auto');
            print(self.fighand,filename,'-dpng',sprintf('-r%d',resolution));
            set(self.fighand,'paperpositionmode',orgppm);
            if self.exportwithannulus
                im = imread(filename);
                ann = self.addannulus(im);
                imwrite(im,filename,'PNG','alpha',ann);
            end
            if self.grayscale
                [im,x,ann] = imread(filename);
                im = uint8(mean(im,3));
                imwrite(im,filename,'PNG','alpha',ann);
            end
        end

        function updateface(self)
        % Update the display of the face based on current camera and
        % lighting properties. Mainly for internal use.
        % updateface()
            % This little nugget prevents updateface from self-calling.
            if isselfcall
                return
            end
            if ~all(ishandle([self.fighand self.axhand self.mhand]))
                % Re-open figure if you've closed the figure or mucked up
                % the axis
                self.openfigure;
            end
            % Ensure visibility is set correctly
            if ~strcmp(get(self.fighand,'visible'),self.visible)
                set(self.fighand,'visible',self.visible);
            end
            % Reset shape/texture data
            if ~all(get(self.mhand,'facevertexcdata')==self.tex)
                set(self.mhand,'facevertexcdata',self.tex);
            end
            if ~all(get(self.mhand,'vertices')==self.shape)
                set(self.mhand,'vertices',self.shape);
            end
            if ~all(get(self.fighand,'color')==self.bgcolor)
                set(self.fighand,'color',self.bgcolor);
            end
            % Set head view
            if ~all(get(self.axhand,'view') == ...
                    [self.azioffset + self.azimuth self.elevation])
                set(self.axhand,'view',...
                    [self.azioffset+self.azimuth self.elevation]);
            end
            % Camera
            if self.centercamoneyes
                camtarget(self.axhand,self.eyesxyz);
            else
                camtarget(self.axhand,self.defxyz);
            end
            camva(self.axhand,self.zoom);
            % ambient lights
            for s = {'ambientstrength','diffusestrength','specularstrength'}
                if ~get(self.mhand,s{1}) == self.(s{1})
                    set(self.mhand,s{1},self.(s{1}));
                end
            end
            % Lighting (maybe none?)
            if isempty(self.lights)
                drawnow;
                return
            end
            % update each light source
            for lind = 1:length(self.lights)
                % Create light source if it isn't already initialised
                if isnan(self.lights(lind).hand)
                    % initialise with dummy azi/el settings (gets updated
                    % later)
                    self.lights(lind).hand = camlight(0,0);
                    set(self.lights(lind).hand,'Color',...
                        self.lights(lind).color,...
                        'Style',self.lights(lind).style);
                end
                % Do we need to comp for head view?
                if self.lights(lind).centeroncam
                    camlight(self.lights(lind).hand,...
                        self.lights(lind).azimuth,...
                        self.lights(lind).elevation);
                else
                    l_azimuth = self.lights(lind).azimuth-self.azimuth;
                    l_elevation = self.lights(lind).elevation- ...
                        self.elevation;
                    camlight(self.lights(lind).hand,l_azimuth,l_elevation);
                end
            end
            drawnow;
            % Update image version
            self.image = self.exportasmatrix;
            if ~isempty(self.imagelims)
                self.image = self.image(...
                    self.imagelims(1):self.imagelims(2),...
                    self.imagelims(3):self.imagelims(4),:);
            end
            % add imagealpha field
            if ~isempty(self.alpha)
                imsize = size(self.image);
                % set background color to self.bgcolor (may be slightly off
                % from 'true' bgcolor from exportasmatrix)
                % ensure double but in uint8 range
                if any(self.bgcolor > 1)
                    bg = double(self.bgcolor) / 255;
                else
                    bg = self.bgcolor;
                end
                if self.grayscale
                    bg = ones(imsize) .* mean(self.bgcolor);
                    alpha = self.alpha;
                else
                    bg = ones(imsize) .* repmat(reshape(self.bgcolor,...
                        [1 1 3]),[imsize(1:2) 1]);
                    % ensure alpha is 3d
                    alpha = repmat(self.alpha,[1 1 imsize(3)]);
                end
                im = double(self.image) / 255;
                self.imagealpha = (im.*alpha) + (bg.*(1-alpha));
                self.imagealpha = im2uint8(self.imagealpha);
            end
        end
    end

    methods (Access = private)

        function openfigure(self)
        % Display a figure. Moved from constructor method to support
        % re-opening a FaceFigure instance later. Called on construction
        % and on updateface if self.fighand goes AWOL. So the figure is
        % only re-opened if you do something with the instance that
        % requires the figure being present.
        % openfigure(self)
            figure(self.fighand);
            set(self.fighand,'Renderer','opengl','name',...
                'FaceFigure','inverthardcopy','off','units','pixels',...
                'visible',self.visible);
            clf;
            fig_pos = get(self.fighand,'Position');
            fig_pos(3:4) = self.figsize;
            set(self.fighand,'Position',fig_pos,'paperpositionmode',...
                'auto');
            set(self.fighand, 'ResizeFcn',@resizeCallback);
            self.axhand = gca;
            % Display mesh
            % nb we now assume that self.shape has already been reordered
            % to (:,[ 1 3 2]).
            self.mhand = trimesh(self.tl,self.shape(:,1),...
                self.shape(:,2),...
                self.shape(:,3),'EdgeColor','none','FaceVertexCData',...
                self.tex,'FaceColor','interp','FaceLighting','phong',...
                'ambientstrength',self.ambientstrength,...
                'diffusestrength',self.diffusestrength,...
                'specularstrength',self.specularstrength,...
                'specularexponent',self.specularexponent);
            % Make display settings sensible
            set(self.axhand,'DataAspectRatio',[1 1 1],...
                'PlotBoxAspectRatio',[1 1 1],'Units','pixels',...
                'GridLineStyle', 'none','Position', [0 0 self.figsize],...
                'Visible', 'off','box','off','Projection','perspective'); 
            set(self.fighand, 'Color', self.bgcolor);
            self.removeposterior(self.trimfaceprop);
            % Store segmentation image and original vertices
            % (but not if it's already been initialised with default)
            if isempty(self.texorg)
                self.texorg = get(self.mhand,'facevertexcdata');
            end
            % Extract eye centre of mass
            if isempty(self.eyesxyz)
                verts = get(self.mhand,'vertices');
                eyeverts = verts(self.segbin(:,2)>0,:);
                self.eyesxyz = mean(eyeverts,1);
            end
            % Store original camtarget to allow dynamic switching between
            % camera modes
            if isempty(self.defxyz)
                self.defxyz = camtarget(self.axhand);
            end
            % Get auto zoom setting
            currzoom = camva(self.axhand);
            % Make manual and possibly zoom
            if isempty(self.zoom)
                self.zoom = currzoom;
            end
            % Now make fig tex grayscale too
            if self.grayscale
                self.changetexture('grayscale');
            end
            % Rotate to default
            self.turnface(self.azimuth,self.elevation);
            % Apply lights if necessary (but not for non-structs to support
            % disabling light completely with e.g. [])
            if isstruct(self.lights) 
                if isempty(self.lights)
                    % add default light
                    self.addlight(self.lightsdefault.azimuth,...
                        self.lightsdefault.elevation)
                else
                    % need to reintroduce lights from old figure (probably
                    % we are opening a new figure window for the instance)
                    oldlights = self.lights;
                    self.lights(:) = [];
                    for l = 1:length(oldlights)
                        cl = oldlights(l);
                        self.addlight(cl.azimuth,cl.elevation,...
                            'style',cl.style,'color',cl.color,...
                            'centeroncam',cl.centeroncam);
                    end
                end
            end
            % setup image alpha
            if self.alphamask
                [self.image,self.alpha] = self.exportasmatrix;
            else
                self.image = self.exportasmatrix;
            end
            % setup imagelims for trimming the image to size
            if self.imageautotrim
                if self.alphamask
                    % use alpha to mask
                    [self.alpha self.imagelims] = imtrim(self.alpha);
                    % this is also in updateface but we need it here to
                    % make sure scrambling gets applied to the cropped
                    % image below
                    self.image = self.image(...
                        self.imagelims(1):self.imagelims(2),...
                        self.imagelims(3):self.imagelims(4),:);
                else
                    % use the image itself
                    [self.image, self.imagelims] = imtrim(self.image);
                end
            end
            self.updateface;
            % phase scramble requires odd-sized image
            im2scr = self.(self.scrambletarget);
            imsz = size(im2scr);
            % pad an extra row or column as needed
            if ~isodd(imsz(1))
                im2scr(end+1,:,:) = im2scr(end,:,:);
            end
            if ~isodd(imsz(2))
                im2scr(:,end+1,:) = im2scr(:,end,:);
            end
            self.scramble = im_phasescramble(im2scr,1);
            % restore to original size
            self.scramble = self.scramble(1:imsz(1),1:imsz(2),:);
        end

        function annulus = addannulus(self,im)
        % add a gaussian annulus to an image.
        % annulus = addannulus(im)
            [yr xr ndim] = size(im);
            xmid = round(xr/2);
            ymid = round(yr/2);

            % Scale annulus size by zoom
            zoomscale = self.zoomfactor/self.zoom;

            % convert scales to pixels
            fh = round(yr * (self.fh*zoomscale));
            fw = round(yr * (self.fw*zoomscale));
            xoffset = round(xr * (self.xoffset*zoomscale));

            % 0 if center cam on eyes, otherwise scale
            if self.centercamoneyes
                yoffset = round(yr * (self.yoffset*zoomscale));
            else
                yoffset = round(yr * (self.yoffsetalt*zoomscale));
            end

            % make a circle
            c = makecircle(1000);
            % scale by face dims
            e = imresize(c,[fh fw]);
            % cut off top
            cuttop = round(fh * (self.cuttop * zoomscale));
            e(1:cuttop,:) = 0;
            es = size(e);
            % Insert in annulus
            annulus = zeros(yr,xr);
            anny = yoffset+(ymid-floor(es(1)/2):yr);
            anny = anny(1:es(1));
            annx = xoffset+(xmid-floor(es(2)/2):xr);
            annx = annx(1:es(2));
            annulus(anny,annx) = e;
            % Smooth
            fsize = round([yr xr] * .05);
            si = mean(fsize) / 3;
            annulus = conv2(annulus,fspecial('gaussian',fsize,si),'same');
        end

    end
end

%% -----------------------------------------CALLBACK--------
% for resizing figure windows with meshes (from display_face)
function resizeCallback (obj, eventdata)
    
    fig = gcbf;
    fig_pos = get(fig, 'Position');

    axis = findobj(get(fig, 'Children'), 'Tag', 'Axis.Head');
    set(axis, 'Position', [ 0 0 fig_pos(3) fig_pos(4) ]);
end
    
