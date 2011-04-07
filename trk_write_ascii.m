function trk_write_ascii(tracks,savePath)
%TRK_WRITE_ASCII - Save a tract group in an ascii format for use in R
%
% Syntax: trk_write_ascii(tracks,savePath)
%
% Inputs:
%    tracks   - Track data in matrix form with 1 scalar attached (assumed to be FA)
%               [nPts x 4 x nTracks]
%    savePath - Path and name of desired output file.
%
% Output files:
%    Saves an ASCII data table with variables: Streamline, Point, x, y, z, FA
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
%    trk_write_ascii(trk_restruc(tracks_sc), 'streamlines.txt')
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: http://github.com/johncolby/along-tract-stats/wiki/single-subject

% Author: John Colby (johncolby@ucla.edu)
% UCLA Developmental Cognitive Neuroimaging Group (Sowell Lab)
% Oct 2010

if nargin < 2 || isempty(savePath); savePath = fullfile(pwd, 'tracks.txt'); end
if nargin < 1 || isempty(tracks); error('Must supply tracks to write.'); end

fid = fopen(savePath, 'w');

nPts    = size(tracks, 1);
nTracks = size(tracks, 3);

outmat = [permute(repmat(1:nTracks,nPts,1), [1 3 2]) repmat((1:nPts)',[1 1 nTracks]) tracks];

fprintf(fid, 'Streamline\tPoint\tx\ty\tz\tFA\n');
fprintf(fid, '%i\t%i\t%f\t%f\t%f\t%0.4f\n', permute(outmat,[2 1 3]));

fclose(fid);
