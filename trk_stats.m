function [meanInt,stdInt,nVox] = trk_stats(header,tracks,volume,method)
%TRK_STATS - Compute intensity statistics about a volume within a track group
%Use this function to extract the mean/etc. of a scalar MRI volume within a
%select track group of interest. For instance, you could calculate the mean
%FA within the left corticospinal tract.
%
% Syntax: [meanInt,stdInt,nVox] = trk_stats(header,tracks,volume,method)
%
% Inputs:
%    header - Header information from .trk file [struc]
%    tracks - Track data struc array [1 x nTracks]
%    volume - Scalar MRI volume
%    method - Method to use for extracting intensities from the MRI volume
%      quick   - Extract intensities directly. This uses linear indexing to  
%                avoid a for loop and is ~ 5x faster than the higher-level
%                nearest neighbor version below.
%      nearest - Nearest neighbor interpolation. Should be equivalent to direct 
%                lookup.
%      cubic   - Cubic interpolation instead of direct voxel lookup.
%
% Outputs:
%    meanInt - Mean image volume intensity (e.g. FA) across an entire tract group
%    stdInt  - Standard deviation of these intensities
%    nVox    - Number of unique voxels traversed by the track group
%
% Example:
%    [header tracks]       = trk_read(trkPath);
%    volume                = read_avw(volPath);
%    [meanInt stdInt nVox] = trk_stats(header, tracks, volume, 'quick')
%
% Other m-files required: trk_restruc
% Subfunctions: none
% MAT-files required: none
%
% See also: TRK_READ, READ_AVW

% Author: John Colby (johncolby@ucla.edu)
% UCLA Developmental Cognitive Neuroimaging Group (Sowell Lab)
% Mar 2010

% Input argument defaults
if nargin < 4 || isempty(method), method = 'quick'; end

% Reformat tracks
vert       = cat(1, tracks.matrix);
tracks_mat = trk_restruc(tracks);
vox        = tracks_mat ./ repmat(repmat(header.voxel_size, size(tracks_mat,1), 1), [1 1 size(tracks_mat,3)]);
vox_rnd    = ceil(vox);

% Extract intensities from MRI volume
if strcmp(method, 'quick')
    inds = sub2ind(header.dim, vox_rnd(:,1,:), vox_rnd(:,2,:), vox_rnd(:,3,:));
    ints = volume(inds(~isnan(inds)));
    meanInt = mean(ints(ints>0));
    stdInt  = std(ints(ints>0));
    nVox    = length(unique(ceil(vert ./ repmat(header.voxel_size, size(vert,1),1)), 'rows'));
else
    vox(isnan(vox)) = 0;
    [x, y, z] = ndgrid(0.5:(header.dim(1)-0.5), 0.5:(header.dim(2)-0.5), 0.5:(header.dim(3)-0.5));
    ints = interpn(x,y,z,volume, vox(:,1,:), vox(:,2,:), vox(:,3,:), method);
    ints(ints==0) = NaN;
    
    % Calculate descriptive statistics
    meanInt = nanmean(ints(ints>0));
    stdInt  = nanstd(ints(ints>0));
    nVox    = length(unique(ceil(vert ./ repmat(header.voxel_size, size(vert,1),1)), 'rows'));
end
