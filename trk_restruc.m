function [tracks_out,tiePoint] = trk_restruc(tracks_in,tiePoint)
%TRK_RESTRUC - Restrucutres streamline data between matrix and structure array
%forms. If all streamlines do not have the same length (i.e. TRK_INTERP has not
%been run) then the matrix will be padded with NaNs.
%
% Syntax: [tracks_out,tiePoint] = trk_restruc(tracks_in,tiePoint)
%
% Inputs:
%    tracks_in - Tracks in matrix form [nPoints x 3+nScalars x nTracks] or
%                Tracks in structure form [1 x nTracks]
%    tiePoint  - Vertex index closest to the midpoint of the mean tract geometry
%                (only relevant if using constant spacing mode for trk_interp)
%
% Outputs:
%    tracks_out - Track data in the opposite form
%    tiePoint   - Modified tiePoint after flipping
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
% See also: TRK_READ, TRK_INTERP

% Author: John Colby (johncolby@ucla.edu)
% UCLA Developmental Cognitive Neuroimaging Group (Sowell Lab)
% Mar 2010

if nargin < 2, tiePoint = []; end

if isstruct(tracks_in) % Structure --> matrix
    nPoints = [tracks_in.nPoints];
    [maxPoints ind] = max(nPoints);
    tracks_out = nan(maxPoints, size(tracks_in(1).matrix, 2), length(tracks_in));
    if length(unique(nPoints))~=1
        
        % If there is a tiePoint given, center the streamlines on this point
        if isfield(tracks_in, 'tiePoint')
            offset = max([tracks_in.tiePoint]) - tracks_in(ind).tiePoint;
            for i=1:length(tracks_in)
                indStart = 1 + offset + (tracks_in(ind).tiePoint - tracks_in(i).tiePoint);
                tracks_out(indStart:(indStart+tracks_in(i).nPoints-1), 1:size(tracks_in(i).matrix, 2), i) = tracks_in(i).matrix;
            end
            tracks_out(tracks_out==0) = NaN;
            tiePoint = tracks_in(ind).tiePoint + offset;
        
        % If there isn't a tiePoint given, start all streamlines in row 1 and
        % only NaN pad the end
        else
            for i=1:length(tracks_in)
                tracks_out(1:size(tracks_in(i).matrix, 1), 1:size(tracks_in(i).matrix, 2), i) = tracks_in(i).matrix;
            end
            tiePoint = [];
        end
        
    % If the streamlines are all the same length, simple algebra will work    
    else
        tracks_out = cat(3, tracks_in.matrix);
        tiePoint = [];
    end
    
else % Matrix --> structure
    [tracks_out(1:size(tracks_in, 3)).nPoints] = deal([]);
    for iTrk=1:size(tracks_in, 3)
        tracks_out(iTrk).matrix  = tracks_in(:,:,iTrk);
        tracks_out(iTrk).nPoints = size(tracks_in, 1);
    end
    if tiePoint
        [tracks_out.tiePoint] = deal(tiePoint);
    end
end
