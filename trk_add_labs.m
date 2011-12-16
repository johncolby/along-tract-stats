function [header,tracks,info] = trk_add_labs(header,tracks,proto,spacing,method,opts)
%TRK_ADD_LABS - Match fiber vertices to vertices on a prototype fiber
%Instead of prescribing correspondence based on vertex index, alternate methods
%can be used. These correspondence labels are added just like a regular scalar
%to the tract file.
%
% Syntax: [header,tracks,info] = trk_add_labs(header,tracks,proto,spacing,method,opts)
%
% Inputs:
%    header  - .trk file header
%    tracks  - .trk file body (structure form)
%    proto   - Prototype fiber (Ex: trk_restruc(mean(tracks_interp, 3)))
%    spacing - Desired arc length spacing along prototype, in mm. (Default: 4mm)
%    method  - Assignment method [string]
%           DM - Distance map method. A discretized lookup table is generated on
%                a voxel-like grid, with labels indicating the closest vertex on
%                the prototype. If multiple tract groups share the same
%                prototype, this can dramatically speed up future processing.
%                (Maddah 2008, PMID: 18180197)
%           OP - Optimal point match method. Instead of total Euclidean
%                distance, only the component tangent to the prototype fiber is
%                considered. Assignment is optimized via the Hungarian
%                algorithm. (O'Donnell 2009, PMID: 19154790)
%    opts    - Method-specific options. [struc]
%           DM - grid_spacing: Grid spacing for distance/label maps
%                (Default: 1mm)
%              - max_endpts: Maximum number of matches allowed to the prototype
%                endpoints (Default: Inf)
%              - info: (optional) Previous distance map info to be reused
%           OP - tan_thresh: A threshold on the allowable distance for matching
%                in the tangent direction, as a fraction of the arc length
%                spacing along the prototype. (Default: 0.4)
%              - euc_thresh: Since OP only considers the tangent component of
%                distance, it could erroneously match far away vertices if a 
%                normal plane of the prototype crosses the fiber group in
%                multiple places (e.g. curved tracts like the AF). Using a
%                Euclidean distance treshold (as a fraction of prototype fiber
%                length) prevents that. (Default: 0.3)
%
% Outputs:
%    header - Updated header
%    tracks - Updated tracks structure
%    info   - (optional) Method-specifc info [struc]
%          DM - dm: Distance map [3D matrix]
%             - lm: Label map [3D matrix]
%             - grid_spacing: Grid spacing for distance/label maps [num]
%             - max_endpts: Max allowable matches to prototype endpoints [num]
%             - xgrid: x grid edges (not centers) [size(dm,1)+1]
%             - ygrid: y grid edges
%             - zgrid: z grid edges
%          OP - not used
%
% Example: 
%    exDir           = '/path/to/along-tract-stats/example';
%    subDir          = fullfile(exDir, 'subject1');
%    trkPath         = fullfile(subDir, 'CST_L.trk');
%    [header tracks] = trk_read(trkPath);
%    tracks_interp   = trk_interp(tracks, 100);
%    track_mean      = mean(tracks_interp, 3);
%    track_mean_str  = trk_restruc(track_mean);
%
%    opts = struct('grid_spacing', 1,...
%                  'max_endpts', 1);
%    [header_dm tracks_dm info] = trk_add_labs(header, tracks, track_mean_str, 4, 'DM', opts);
%
% Other m-files required: assignmentoptimal.m, join.m
% Subfunctions: doDM, doOP, biggest_segment
% MAT-files required: none
%
% See also: TRK_READ, TRK_INTERP, ASSIGNMENTOPTIMAL, JOIN

% Author: John Colby (johncolby@ucla.edu)
% UCLA Developmental Cognitive Neuroimaging Group (Sowell Lab)
% Oct 2011

if nargin < 6, opts = struct(); end
if nargin < 5 || isempty(method), error('Must specify a method.'); end
if isempty(spacing), spacing = 4; end
info = [];

% Linear interpolation of prototype to arc length coordinates
segs = sqrt(sum((proto.matrix(2:end,1:3) - proto.matrix(1:(end-1),1:3)).^2, 2));
dist = [0; cumsum(segs)];
proto.matrix = interp1(dist, proto.matrix, min(dist):spacing:max(dist), 'linear');

% Obtain labels matching fiber vertices to prototype vertices
switch method
    
    % Distance map (Maddah 2008)
    case 'DM'
        if ~isfield(opts, 'grid_spacing') || isempty(opts.grid_spacing), opts.grid_spacing = 1; end
        if ~isfield(opts, 'max_endpts') || isempty(opts.max_endpts), opts.max_endpts = Inf; end
        if ~isfield(opts, 'info'), opts.info = []; end
        
        [tracks info] = doDM(tracks, proto, opts.grid_spacing, opts.max_endpts, opts.info);
    
    % Optimal point match (O'Donnell 2009)
    case 'OP'
        if ~isfield(opts, 'tan_thresh') || isempty(opts.tan_thresh), opts.tan_thresh = 0.4; end
        if ~isfield(opts, 'euc_thresh') || isempty(opts.euc_thresh), opts.euc_thresh = 0.3; end
        
        opts.tan_thresh = opts.tan_thresh * spacing;
        opts.euc_thresh = opts.euc_thresh * trk_length(proto);
        
        tracks = doOP(tracks, proto, opts.tan_thresh, opts.euc_thresh);
end

% Update header
name = join('_', 'Label', method);
n_scalars_old    = header.n_scalars;
header.n_scalars = n_scalars_old + 1;
header.scalar_name(n_scalars_old + 1,1:size(name,2)) = name;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Subfunctions

% Distance map method
function [tracks,info] = doDM(tracks,proto,grid_spacing,max_endpts,info)

if(isempty(info))
    % Define grid for distance map
    tracks_mat = trk_restruc(tracks);
    
    x = tracks_mat(:,1,:);
    y = tracks_mat(:,2,:);
    z = tracks_mat(:,3,:);
    
    xgrid = floor(min(x(:))):grid_spacing:ceil(max(x(:)));
    ygrid = floor(min(y(:))):grid_spacing:ceil(max(y(:)));
    zgrid = floor(min(z(:))):grid_spacing:ceil(max(z(:)));
    
    [X,Y,Z] = meshgrid(xgrid(1:(end-1)), ygrid(1:(end-1)), zgrid(1:(end-1)));
    centers = [X(:) Y(:) Z(:)] + grid_spacing/2;
    
    % Calculate distance and label maps
    % Alternatively, it is much faster to loop over prototype vertices, if you
    % have enough memory to operate on all voxels at once
    [dists, labs] = deal(zeros(length(centers), 1));
    for ivox=1:length(centers)
        [dists(ivox), labs(ivox)] = min(sqrt(sum(bsxfun(@minus, centers(ivox,:), proto.matrix(:,1:3)).^2, 2)));
    end
    
    % Save assignment info for future reuse
    info.dm           = reshape(dists, size(X));
    info.lm           = reshape(labs, size(X));
    info.grid_spacing = grid_spacing;
    info.max_endpts   = max_endpts;
    info.xgrid        = xgrid;
    info.ygrid        = ygrid;
    info.zgrid        = zgrid;
end

% Assign vertices of other fibers to vertices on the prototype
for iTrk=1:length(tracks)
    inds = ceil(bsxfun(@minus, tracks(iTrk).matrix(:,1:3), [min(info.xgrid) min(info.ygrid) min(info.zgrid)]) / info.grid_spacing);
    labs = info.lm(sub2ind(size(info.lm), inds(:,2), inds(:,1), inds(:,3)));
    
    % Limit the number of fiber points that match the prototype endpoints
    if length(labs(labs==labs(1))) > max_endpts
        labs(labs==labs(1))   = [repmat(NaN, length(labs(labs==labs(1)))-max_endpts, 1);...
                                 repmat(labs(1), max_endpts, 1)];
    end
    if length(labs(labs==labs(end))) > max_endpts
    labs(labs==labs(end)) = [repmat(labs(end), max_endpts, 1);...
                             repmat(NaN, length(labs(labs==labs(end)))-max_endpts, 1)];
    end
    
    % Trim labs to only the largest monotonic segment
    labs = biggest_segment(labs);
    
    tracks(iTrk).matrix = [tracks(iTrk).matrix labs];
end

% Optimal point match method
function tracks = doOP(tracks,proto,tan_thresh,euc_thresh)

for iTrk=1:length(tracks)
    % Calculate displacement vectors between all vertex pairs
    [X Y] = meshgrid(1:size(tracks(iTrk).matrix,1), 1:size(proto.matrix,1));
    v = proto.matrix(Y(:),1:3) - tracks(iTrk).matrix(X(:),1:3);
    
    % Calculate Euclidean distances
    D_euc = sqrt(sum(v.^2, 2));
    distmat_euc = reshape(D_euc, size(X));
    
    % Calculate tangent tensors along prototype
    tangents = diff(proto.matrix(:,1:3));
    tangents = bsxfun(@rdivide, tangents, sqrt(sum(tangents.^2, 2))); %normalize to unit length
    for iseg=1:length(tangents)
        M(:,:,iseg) = tangents(iseg,:)' * tangents(iseg,:);
    end
    M = cat(3, M, M(:,:,end));
    
    % Calculate distance components that are along the tangents
    % Egn. 1 in O'Donnell 2009
    D = zeros(length(v),1);
    for iv=1:length(v)
        D(iv) = sqrt(v(iv,:) * M(:,:,Y(iv)) * v(iv,:)');
    end
    distmat = reshape(D, size(X));
    
    % Tangential distance threshold
    distmat(distmat > tan_thresh) = Inf;
    
    % Euclidean distance threshold
    distmat(distmat_euc > euc_thresh) = Inf;
    
    % Assign vertices with Hungarian algorithm
    labs = assignmentoptimal(distmat');
    
    % Linear interpolation followed by rounding for fiber vertices that didn't
    % match the prototype
    labs(labs==0) = NaN;
    labs = round(interp1(find(~isnan(labs)), labs(~isnan(labs)), 1:length(labs), 'linear'))';
    
    % Trim labs to only the largest monotonic segment
    labs = biggest_segment(labs);
    
    tracks(iTrk).matrix = [tracks(iTrk).matrix labs];
end

% Trim labels to only the largest monotonically increasing or decreasing segment
function labs = biggest_segment(labs)

dlabs = diff(labs);
if mean(dlabs(~isnan(dlabs))) < 0
    dlabs = -dlabs;
end
segs = [0; cumsum(dlabs < 0)];
labs(segs ~= mode(segs)) = NaN;