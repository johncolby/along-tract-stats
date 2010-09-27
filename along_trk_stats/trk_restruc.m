function tracks_out = trk_restruc(tracks_in)
%TRK_RESTRUC - Restrucutres track data between matrix and strucutre array forms
%If all tracks do not have the same length (i.e. TRK_INTERP has not been run)
%then the matrix will be padded with NaNs.
%
% Syntax: tracks_out = trk_restack(tracks_in)
%
% Inputs:
%    tracks_in - Tracks in matrix form [nPoints x 3+nScalars x nTracks] or
%                Tracks in structure form [1 x nTracks]
%
% Outputs:
%    tracks_out - Track data in the opposite form
%
% Example:
%    [header tracks]   = trk_read(filePath)
%    tracks_interp     = trk_interp(tracks, 100);
%    tracks_interp_str = trk_restruc(tracks_interp);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: TRK_READ, TRK_INTERP

% Author: John Colby (johncolby@ucla.edu)
% UCLA Developmental Cognitive Neuroimaging Group (Sowell Lab)
% Mar 2010 $Rev$ $Date$

if isstruct(tracks_in) % Structure --> matrix
    nPoints = cat(tracks_in.nPoints);
    if length(unique(nPoints))~=1
        warning('Tracks don''t have same length; padded with NaNs.')
        tracks_out = nan(max(nPoints), size(tracks_in(1).matrix, 2), length(tracks_in));
        for i=1:length(tracks_in)
            tracks_out(1:size(tracks_in(i).matrix, 1), 1:size(tracks_in(i).matrix, 2), i) = tracks_in(i).matrix;
        end
    else
        tracks_out = cat(3, tracks_in.matrix);
    end
else % Matrix --> structure
    [tracks_out(1:size(tracks_in, 3)).nPoints] = deal([]);
    for iTrk=1:size(tracks_in, 3)
        mat                       = tracks_in(:,:,iTrk);
        mat(all(isnan(mat), 2),:) = [];
        tracks_out(iTrk).matrix   = mat;
        tracks_out(iTrk).nPoints  = size(mat, 1);
    end
end