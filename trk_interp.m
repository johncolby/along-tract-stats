function [tracks_interp,trk_mean_length] = trk_interp(tracks,nPoints_new,spacing)
%TRK_INTERP - Interpolate tracks with cubic B-splines
%Each interpolated track polyline will have the same number of vertices. Note
%however that this does not give uniform spatial sampling, so shorter tracks
%will end up with denser sampling than longer tracks. May be useful to groom
%your tracks first. Can be multithreaded if Parallel Computing Toolbox is
%installed.
%
% Syntax: [tracks_interp,trk_mean_length] = trk_interp(tracks,nPoints_new,spacing)
%
% Inputs:
%    tracks      - Struc array output of TRK_READ [1 x nTracks]
%    nPoints_new - Number of vertices for each streamlines (spacing 
%                  between vertices will vary between streamlines)
%    spacing     - Spacing between each vertex (# of vertices will vary between
%                  streamlines). Note: Only supply nPoints_new *OR* spacing!
%
% Outputs:
%    tracks_interp   - Interpolated tracks [nPoints_new x 3 x nTracks]
%    trk_mean_length - The length of the mean tract geometry if using the
%                      spacing parameter, above. Useful to normalize track
%                      lengths between subjects.
%
% Example:
%    trkPath                 = fullfile(exDir, 'cst.trk');
%    [header tracks]         = trk_read(trkPath);
%    tracks_interp           = trk_interp(tracks, 100);
%    tracks_interp           = trk_flip(header, tracks_interp, [97 110 4]);
%    tracks_interp_str       = trk_restruc(tracks_interp);
%    [header_sc tracks_sc]   = trk_add_sc(header, tracks_interp_str, volume, 'FA');
%    [scalar_mean scalar_sd] = trk_mean_sc(header_sc, tracks_sc);
%
% Other m-files required: Spline Toolbox
% Subfunctions: none
% MAT-files required: none
%
% See also: TRK_READ, SPLINE, PARFOR

% Author: John Colby (johncolby@ucla.edu)
% UCLA Developmental Cognitive Neuroimaging Group (Sowell Lab)
% Apr 2010

if nargin<3, spacing = []; end
if nargin<2 || isempty(nPoints_new), nPoints_new = 100; end

tracks_interp   = zeros(nPoints_new, 3, length(tracks));
trk_mean_length = [];
pp = repmat({[]},length(tracks));

% Interpolate streamlines so that each has the same number of vertices, spread
% evenly along each streamline (i.e. vertex spacing will vary between streamlines)
parfor iTrk=1:length(tracks)
    tracks_tmp = tracks(iTrk);
    
    % Determine streamline segment lengths
    segs = sqrt(sum((tracks_tmp.matrix(2:end,1:3) - tracks_tmp.matrix(1:(end-1),1:3)).^2, 2));
    dist = [0; cumsum(segs)];
    
    % Fit spline
    pp{iTrk} = spline(dist, tracks_tmp.matrix');
    
    % Resample streamline along the spline
    tracks_interp(:,:,iTrk) = ppval(pp{iTrk}, linspace(0, max(dist), nPoints_new))';
end

if ~isempty(spacing)
    % Calculate streamline lengths
    lengths = trk_length(tracks_interp);
    
    % Determine the mean tract geometry and grab the middle vertex
    track_mean      = mean(tracks_interp, 3);
    trk_mean_length = trk_length(track_mean);
    middle          = track_mean(round(length(track_mean)/2),:);
    
    % Interpolate streamlines again, but this time sample with constant vertex
    % spacing for all streamlines. This means that the longer streamlines will now
    % have more vertices.
    tracks_interp = repmat(struct('nPoints', 0, 'matrix', [], 'tiePoint', 0), 1, length(tracks));
    parfor iTrk=1:length(tracks)
        tracks_interp(iTrk).matrix  = ppval(pp{iTrk}, 0:spacing:lengths(iTrk))';
        tracks_interp(iTrk).nPoints = size(tracks_interp(iTrk).matrix, 1);
        
        % Also determine which vertex is the "tie down" point by finding the one
        % closest to the middle point of the mean tract geometry
        dists = sqrt(sum(bsxfun(@minus, tracks_interp(iTrk).matrix, middle).^2,2));
        [~, ind] = min(dists);
        tracks_interp(iTrk).tiePoint = ind;
    end
end
