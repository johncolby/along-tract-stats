function [scalar_mean,scalar_sd] = trk_mean_sc(header,tracks)
%TRK_MEAN_SC - Calculate the mean scalar along a track
%Returns the mean and SD of a scalar volume (e.g. FA map) *along* a track.
%Rather than collapsing across the whole track, as in TrackVis or TRK_STATS,
%this function returns vectors corresponding to the different vertices along the
%whole track. This will allow you to localize differences within a track.
%
% Syntax: [scalar_mean,scalar_sd] = trk_mean_sc(header,tracks)
%
% Inputs:
%    header - Header information from .trk file [struc]
%    tracks - Track data struc array [1 x nTracks]
%
% Outputs:
%    scalar_mean - Mean of the scalar at each track point [nPoints x nScalars]
%    scalar_sd   - Standard deviation of the scalar at each track point
%                  [nPoints x nScalars]
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
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: TRK_READ, READ_AVW, TRK_INTERP, TRK_RESTRUC, TRK_ADD_SC

% Author: John Colby (johncolby@ucla.edu)
% UCLA Developmental Cognitive Neuroimaging Group (Sowell Lab)
% Apr 2010

scalars = zeros(tracks(1).nPoints, header.n_count, header.n_scalars);

for i=1:header.n_scalars
    mat_long        = cat(1, tracks.matrix);
    scalars(:,:,i)  = reshape(mat_long(:,4), tracks(1).nPoints, header.n_count, header.n_scalars);
end

scalar_mean = squeeze(nanmean(scalars, 2));
scalar_sd   = squeeze(nanstd(scalars, 0, 2));
