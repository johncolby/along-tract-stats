function [tracks_out pt_start] = trk_flip(header,tracks_in,pt_start,volume,slices)
%TRK_FLIP - Flip the ordering of tracks
%When TrackVis stores .trk files, the ordering of the points are not always
%optimal (e.g. the corpus callosum will have some tracks starting on the left
%and some on the right). TRK_FLIP solves this problem by reordering streamlines
%so that the terminal points nearest to point 'pt_start' will be the starting
%points.
%
% Syntax: [tracks_out,pt_start] = trk_flip(header,tracks_in,pt_start,volume,slices)
%
% Inputs:
%    header    - Header information from .trk file [struc]
%    tracks_in - Tracks in matrix or structure form. Should NOT be padded with
%                NaNs.
%    pt_start  - XYZ voxel coordinates to which streamline start points will be
%                matched. If not given, will determine interactively. [1 x 3]
%    volume    - (optional) Useful if determining pt_start interactively
%    slices    - (optional) Slice planes for 'volume'
%
% Outputs:
%    tracks_out - Output tracks (matrix or struc form, depending on input). Same
%                 vertices as tracks_in, but the ordering of some tracks will
%                 now be reversed.
%    pt_start   - Useful to collect the interactively found pt_starts.
%
% Example:
%    exDir                   = '/path/to/along-tract-stats/example';
%    subDir                  = fullfile(exDir, 'subject1');
%    trkPath                 = fullfile(subDir, 'CST_L.trk');
%    volPath                 = fullfile(subDir, 'dti_fa.nii.gz');
%    volume                  = read_avw(volPath);
%    [header tracks]         = trk_read(trkPath);
%    tracks_interp           = trk_interp(tracks, 100);
%    tracks_interp           = trk_flip(header, tracks_interp, [97 110 4]);
%    tracks_interp_str       = trk_restruc(tracks_interp);
%    [header_sc tracks_sc]   = trk_add_sc(header, tracks_interp_str, volume, 'FA');
%    [scalar_mean scalar_sd] = trk_mean_sc(header_sc, tracks_sc);
%
% Other m-files required: trk_plot, trk_restruc
% Subfunctions: none
% MAT-files required: none
%
% See also: TRK_READ, TRK_INTERP
%           http://www.colbyimaging.com/wiki/neuroimaging/along-tract-stats#interactive_trk_flip

% Author: John Colby (johncolby@ucla.edu)
% UCLA Developmental Cognitive Neuroimaging Group (Sowell Lab)
% Apr 2010

if nargin < 5, slices = []; end
if nargin < 4, volume = []; end

% If no pt_start given, determine interactively
if nargin < 3 || isempty(pt_start) || any(isnan(pt_start))
    if isnumeric(tracks_in)
        tracks_in_str = trk_restruc(tracks_in);
    else
        tracks_in_str = tracks_in;
    end
    fh = figure;
    trk_plot(header, tracks_in_str, volume, slices);
    dcm_obj = datacursormode(fh);
    datacursormode(fh, 'on')
    set(fh,'DeleteFcn','global c_info, c_info = getCursorInfo(datacursormode(gcbo));')
    waitfor(fh)
    global c_info
    pt_start = c_info.Position;
end

tracks_out = tracks_in;

% Fast algebra if streamlines are all the same length
if isnumeric(tracks_in)
    if any(isnan(tracks_in(:)))
        error('If you are going to deal with streamlines padded with NaNs (i.e. different lengths), they should be flipped FIRST.')
    end
    % Determine if the first or last track point is closer to 'pt_start'
    if length(size(tracks_in))==2  % Only 1 streamliine
        point_1   = sqrt(sum(bsxfun(@minus, tracks_in(1,:,:), pt_start).^2, 2));
        point_end = sqrt(sum(bsxfun(@minus, tracks_in(end,:,:), pt_start).^2, 2));
    else
        point_1   = sqrt(sum(squeeze(bsxfun(@minus, tracks_in(1,:,:), pt_start))'.^2, 2));
        point_end = sqrt(sum(squeeze(bsxfun(@minus, tracks_in(end,:,:), pt_start))'.^2, 2));
    end
    
    % Flip the tracks whose first points are not closest to 'pt_start'
    ind                 = point_end < point_1;
    tracks_out(:,:,ind) = tracks_in(fliplr(1:end),:,ind);

% Otherwise, loop through one by one
else
    if any(isnan(cat(1,tracks_in.matrix)))
        error('If you are going to deal with streamlines padded with NaNs (i.e. different lengths), they should be flipped FIRST.')
    end
    if size(tracks_in(1).matrix, 2) > 3
       error('Streamlines should be flipped before scalars are attached.') 
    end
    for iTrk=1:length(tracks_in)
        % Determine if the first or last track point is closer to 'pt_start'
        point_1   = sqrt(sum((tracks_in(iTrk).matrix(1,:) - pt_start).^2));
        point_end = sqrt(sum((tracks_in(iTrk).matrix(end,:) - pt_start).^2));
        
        % Flip the tracks whose first points are not closest to 'pt_start'
        if point_end < point_1;
            tracks_out(iTrk).matrix   = flipud(tracks_in(iTrk).matrix);
            if isfield(tracks_out, 'tiePoint')
                tracks_out(iTrk).tiePoint = tracks_out(iTrk).nPoints - (tracks_out(iTrk).tiePoint-1);
            end
        end
    end
end
