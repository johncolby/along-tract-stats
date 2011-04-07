function [meanInt,stdInt,nVox] = trk_stats_quick(header,tracks,volume,method)
%TRK_STATS_QUICK - Compute intensity statistics about a volume within a track group
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
%                lookup. Only implemented to test code for bugs.
%      cubic   - Cubic interpolation instead of direct voxel lookup.
%
% Outputs:
%    meanInt - Mean image volume intensity (e.g. FA) along a track group
%    stdInt  - Standard deviation of these intensities
%    nVox    - Number of unique voxels traversed by the track group
%
% Example:
%    [header tracks]       = trk_read(trkPath);
%    volume                = read_avw(volPath);
%    [meanInt stdInt nVox] = trk_stats(header, tracks, volume, 'quick')
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: TRK_READ, READ_AVW

% Author: John Colby (johncolby@ucla.edu)
% UCLA Developmental Cognitive Neuroimaging Group (Sowell Lab)
% Mar 2010

% Input argument defaults
if nargin < 4 || isempty(method), method = 'quick'; end

% Translate continuous vertex coordinates into discrete voxel coordinates
vert         = cat(1, tracks.matrix);
vert_vox     = vert ./ repmat(header.voxel_size, size(vert,1),1);
vert_vox_rnd = ceil(vert_vox);

% Extract intensities from MRI volume
if strcmp(method, 'quick')
    inds = sub2ind(header.dim, vert_vox_rnd(:,1), vert_vox_rnd(:,2), vert_vox_rnd(:,3));
    ints = volume(inds);
else
    [x, y, z] = ndgrid(0.5:(header.dim(1)-0.5), 0.5:(header.dim(2)-0.5), 0.5:(header.dim(3)-0.5));
    ints      = interpn(x, y, z, volume, vert_vox(:,1), vert_vox(:,2), vert_vox(:,3), method);
end

% Calculate descriptive statistics
meanInt = mean(ints(ints>0));
stdInt  = std(ints(ints>0));
nVox    = length(unique(vert_vox_rnd, 'rows'));
