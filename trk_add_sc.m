function [header,tracks] = trk_add_sc(header,tracks,volume,name)
%TRK_ADD_SC - Attaches a scalar value to each vertex in a .trk track group
%For example, this function can look in an FA volume, and attach the
%corresponding voxel FA value to each streamline vertex.
%
% Syntax: [header,tracks] = trk_add_sc(header,tracks,volume,name)
%
% Inputs:
%    header - Header information from .trk file [struc]
%    tracks - Track data struc array [1 x nTracks]
%    volume - Scalar MRI volume
%    name   - Description of the scalar to add to the header (e.g. 'FA')
%
% Outputs:
%    header - Updated header
%    tracks - Updated tracks structure
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
% MAT-files required: TRK_READ, READ_AVW
%
% See also: 

% Author: John Colby (johncolby@ucla.edu)
% UCLA Developmental Cognitive Neuroimaging Group (Sowell Lab)
% Apr 2010

% Loop over # of tracks (slow...any faster way?)
for iTrk=1:length(tracks)
    % Translate continuous vertex coordinates into discrete voxel coordinates
    vox = ceil(tracks(iTrk).matrix(:,1:3) ./ repmat(header.voxel_size, tracks(iTrk).nPoints,1));
    
    % Index into volume to extract scalar values
    inds                = sub2ind(header.dim, vox(:,1), vox(:,2), vox(:,3));
    scalars             = volume(inds);
    tracks(iTrk).matrix = [tracks(iTrk).matrix, scalars];
end

% Update header
n_scalars_old    = header.n_scalars;
header.n_scalars = n_scalars_old + 1;
header.scalar_name(n_scalars_old + 1,1:size(name,2)) = name;