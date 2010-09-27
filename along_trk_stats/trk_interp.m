function tracks_interp = trk_interp(tracks,nPoints_new)
%TRK_INTERP - Interpolate tracks with cubic B-splines
%Each interpolated track polyline will have the same number of vertices. Note
%however that this does not give uniform spatial sampling, so shorter tracks
%will end up with denser sampling than longer tracks. May be useful to groom
%your tracks first.
%
% Syntax: tracks_interp = trk_interp(tracks,nPoints_new)
%
% Inputs:
%    tracks      - Struc array output of TRK_READ [1 x nTracks]
%    nPoints_new - Number of interpolation points for each track
%
% Outputs:
%    tracks_interp - Interpolated tracks [nPoints_new x 3 x nTracks]
%
% Example:
%    [header tracks] = trk_read(filePath)
%    tracks_interp   = trk_interp(tracks, 100);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: TRK_READ

% Author: John Colby (johncolby@ucla.edu)
% UCLA Developmental Cognitive Neuroimaging Group (Sowell Lab)
% Apr 2010 $Rev$ $Date$

tracks_interp = zeros(nPoints_new, 3, length(tracks));

parfor iTrk=1:length(tracks)
    tracks_tmp = tracks(iTrk);
    segs = sqrt(sum((tracks_tmp.matrix(2:end,1:3) - tracks_tmp.matrix(1:(end-1),1:3)).^2, 2));
    dist = [0; cumsum(segs)];
    pp   = spline(dist, tracks_tmp.matrix');
    tracks_interp(:,:,iTrk) = ppval(pp, linspace(0, max(dist), nPoints_new))';
end