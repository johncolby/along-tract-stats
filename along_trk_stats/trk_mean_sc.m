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
%    scalar_mean - Mean of the scalar at each track point [nx1]
%    scalar_sd   - Standard deviation of the scalar at each track point [nx1]
%
% Example: 
%    [header tracks]         = read_trk(trkPath);
%    volume                  = read_avw(volPath);
%    tracks_interp           = trk_interp(tracks, 100);
%    tracks_interp_str       = trk_restruc(tracks_interp);
%    [header_sc tracks_sc]   = trk_add_sc(header,tracks_interp_str,volume,'FA')
%    [scalar_mean scalar_sd] = trk_mean_sc(header_sc,tracks_sc);
%    errorbar(scalar_mean, scalar_sd)
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: TRK_READ, READ_AVW, TRK_INTERP, TRK_RESTRUC, TRK_ADD_SC

% Author: John Colby (johncolby@ucla.edu)
% UCLA Developmental Cognitive Neuroimaging Group (Sowell Lab)
% Apr 2010 $Rev$ $Date$

for i=1:header.n_scalars
    mat_long = cat(1, tracks.matrix);
    scalars  = reshape(mat_long(:,4), tracks(1).nPoints, header.n_count);
end

scalar_mean = mean(scalars, 2);
scalar_sd   = std(scalars, 0, 2);