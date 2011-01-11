function trk_write_ascii(tracks,filePath)
%FUNCTION_NAME - One line description of what the function or script performs (H1 line)
%Optional file header info (to give more details about the function than in the H1 line)
%Optional file header info (to give more details about the function than in the H1 line)
%Optional file header info (to give more details about the function than in the H1 line)
%
% Syntax: [output1,output2] = function_name(input1,input2,input3)
%
% Inputs:
%    input1 - Description
%    input2 - Description
%    input3 - Description
%
% Outputs:
%    output1 - Description [1 x n]
%    output2 - Description strc array
%      field1 - Description
%      field2 - Description
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

if nargin < 2 || isempty(filePath); filePath = fullfile(pwd, 'tracks.txt'); end
if nargin < 1 || isempty(tracks); error('Must supply tracks to write.'); end

fid = fopen(filePath, 'w');

nPts    = size(tracks, 1);
nTracks = size(tracks, 3);

outmat = [permute(repmat(1:nTracks,nPts,1), [1 3 2]) repmat((1:nPts)',[1 1 nTracks]) tracks];

fprintf(fid, 'Streamline\tPoint\tx\ty\tz\tFA\n');
fprintf(fid, '%i\t%i\t%f\t%f\t%f\t%0.4f\n', permute(outmat,[2 1 3]));

fclose(fid);
