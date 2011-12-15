function trk_write_ascii(header,tracks,savePath)
%TRK_WRITE_ASCII - Save a tract group in an ascii format for use in R or other
%external programs
%
% Syntax: trk_write_ascii(header,tracks,savePath)
%
% Inputs:
%    tracks   - .trk file body (structure form)
%    savePath - Path and name of desired output file.
%
% Output files:
%    Saves an ASCII data table with variables: Streamline, Point, x, y, z, and
%    (optionally) 1 or more scalars
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
%    trk_write_ascii(tracks_sc, 'streamlines.txt')
%
% Other m-files required: join.m
% Subfunctions: none
% MAT-files required: none
%
% See also: JOIN
%           http://github.com/johncolby/along-tract-stats/wiki/single-subject

% Author: John Colby (johncolby@ucla.edu)
% UCLA Developmental Cognitive Neuroimaging Group (Sowell Lab)
% Oct 2010

if nargin < 2 || isempty(savePath); savePath = fullfile(pwd, 'tracks.txt'); end
if nargin < 1 || isempty(tracks); error('Must supply tracks to write.'); end

fid = fopen(savePath, 'w');

tracks = trk_restruc(tracks);

nPts    = size(tracks, 1);
nTracks = size(tracks, 3);

outmat = [permute(repmat(1:nTracks,nPts,1), [1 3 2]) repmat((1:nPts)',[1 1 nTracks]) tracks];

% Check for scalars
scalar_names = {};
for iscalar = 1:header.n_scalars
    scalar_names{iscalar} = strcat(header.scalar_name(iscalar,:));
end

% Print header row
format = 'Streamline\tPoint\tx\ty\tz';
if ~isempty(scalar_names)
    scalar_names = join('\t', scalar_names);
    format = [format '\t' scalar_names];
end
fprintf(fid, [format '\n']);

% Print along-tract data
format = '%i\t%i\t%f\t%f\t%f';
if ~isempty(scalar_names)
    format = [format '\t' join('\t', repmat({'%f'}, 1, header.n_scalars))];
end
fprintf(fid, [format '\n'], permute(outmat,[2 1 3]));

fclose(fid);
