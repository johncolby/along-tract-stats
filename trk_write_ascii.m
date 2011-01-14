function trk_write_ascii(tracks,savePath)
%TRK_WRITE_ASCII - Save a tract group in an ASCII format for plotting in R
%
% Syntax: trk_write_ascii(tracks,savePath)
%
% Inputs:
%    tracks   - 
%    savePath - 
%
% Output files:
%    Saves an ASCII data table with variables: Streamline, Point, x, y, z, FA
%
% Example: 
%    Line 1 of example
%    Line 2 of example
%    Line 3 of example
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: OTHER_FUNCTION_NAME1,  OTHER_FUNCTION_NAME2

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
