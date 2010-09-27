function [tracks_out pt_start] = trk_flip(header,tracks_in,pt_start,volume)
%TRK_FLIP - Flip the ordering of tracks
%When TrackVis stores .trk files, the ordering of the points are not always
%optimal (e.g. the corpus callosum will have some tracks starting on the left
%and some on the right). TRK_FLIP attempts to help this problem by reordering
%tracks so that the terminal points nearest to point 'pt_start' will be the
%starting points.
%
% Syntax: [tracks_out pt_start] = trk_flip(header,tracks_in,pt_start,volume)
%
% Inputs:
%    header    - Header information from .trk file [struc]
%    tracks_in - Tracks in matrix form. All should have equal number of points.
%    pt_start  - XYZ voxel coordinates to which track start points will be
%                matched. If not given, will determine interactively [1 x 3]
%    volume    - (optional) Useful if determining pt_start interactively
%
% Outputs:
%    tracks_out - Output track matrix. Same vertices as tracks_in, but the
%                 ordering of some tracks will now be reversed.
%    pt_start   - Useful to collect the interactively found pt_starts
%
% Example:
%    Try to get the corpus callosum start points in the left hemisphere
%    tracks_interp_flp = trk_flip(header, tracks_interp, [145 39 28]);
%
% Other m-files required: trk_plot
% Subfunctions: none
% MAT-files required: none
%
% See also: TRK_READ, TRK_INTERP

% Author: John Colby (johncolby@ucla.edu)
% UCLA Developmental Cognitive Neuroimaging Group (Sowell Lab)
% Apr 2010 $Rev$ $Date$

if nargin < 4, volume = []; end

% If no pt_start given, determine interactively
if nargin < 3 || isempty(pt_start) || any(isnan(pt_start))
    fh = figure;
    trk_plot(header, trk_restruc(tracks_in), volume, [48 48 4]);
    dcm_obj = datacursormode(fh);
    datacursormode(fh, 'on')
    set(fh,'DeleteFcn','global c_info, c_info = getCursorInfo(datacursormode(gcbo));')
    waitfor(fh)
    global c_info
    pt_start = c_info.Position;
else
    % Switch from voxels to mm
%     pt_start(2) = header.dim(2) - pt_start(2);   % Temporary tweak until I can
%     pt_start    = pt_start .* header.voxel_size; % figure out the orientation issues
    
end

% Determine if the first or last track point is closer to 'pt_start'
if header.n_count==1
    point_1   = sqrt(sum(bsxfun(@minus, tracks_in(1,:,:), pt_start).^2, 2));
    point_end = sqrt(sum(bsxfun(@minus, tracks_in(end,:,:), pt_start).^2, 2));
else
    point_1   = sqrt(sum(squeeze(bsxfun(@minus, tracks_in(1,:,:), pt_start))'.^2, 2));
    point_end = sqrt(sum(squeeze(bsxfun(@minus, tracks_in(end,:,:), pt_start))'.^2, 2));
end

% Flip the tracks whose first points are not closest to 'pt_start'
ind                 = point_end < point_1;
tracks_out          = tracks_in;
tracks_out(:,:,ind) = tracks_in(fliplr(1:end),:,ind);
