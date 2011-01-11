function [header,tracks] = trk_read(filePath)
%TRK_READ - Load TrackVis .trk files
%
% Syntax: [header,tracks] = trk_read(filePath)
%
% Inputs:
%    filePath - Full path to .trk file [char]
%
% Outputs:
%    header - Header information from .trk file [struc]
%    tracks - Track data struc array [1 x nTracks]
%      nPoints - # of points in each track
%      matrix  - XYZ coordinates and associated scalars [nPoints x 3+nScalars]
%      props   - Properties of the whole tract
%
% Other m-files required: none
% Subfunctions: get_header
% MAT-files required: none
%
% See also: http://www.trackvis.org/docs/?subsect=fileformat

% Author: John Colby (johncolby@ucla.edu)
% UCLA Developmental Cognitive Neuroimaging Group (Sowell Lab)
% Mar 2010

% Parse in header
fid    = fopen(filePath, 'r');
header = get_header(fid);

% Check for byte order
if header.hdr_size~=1000
    fclose(fid);
    fid    = fopen(filePath, 'r', 'b'); % Big endian for old PPCs
    header = get_header(fid);
end

if header.hdr_size~=1000, error('Header length is wrong'), end

% Parse in body
tracks(header.n_count).nPoints = 0;

for iTrk = 1:header.n_count
    tracks(iTrk).nPoints = fread(fid, 1, 'int');
    tracks(iTrk).matrix  = fread(fid, [3+header.n_scalars, tracks(iTrk).nPoints], 'float')';
    if header.n_properties
        tracks(iTrk).props = fread(fid, header.n_properties, 'float');
    end
    
    tracks(iTrk).matrix(:,2) = header.dim(2)*header.voxel_size(2) - tracks(iTrk).matrix(:,2);
    
    % Modify orientation
    %if header.invert_x==1, tracks(iTrk).matrix(:,1) = header.dim(1)*header.voxel_size(1) - tracks(iTrk).matrix(:,1); end
    %if header.invert_y==1, tracks(iTrk).matrix(:,2) = header.dim(2)*header.voxel_size(2) - tracks(iTrk).matrix(:,2); end
    %if header.invert_z==1, tracks(iTrk).matrix(:,3) = header.dim(3)*header.voxel_size(3) - tracks(iTrk).matrix(:,3); end
end

fclose(fid);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function header = get_header(fid)

header.id_string                 = fread(fid, 6, '*char')';
header.dim                       = fread(fid, 3, 'short')';
header.voxel_size                = fread(fid, 3, 'float')';
header.origin                    = fread(fid, 3, 'float')';
header.n_scalars                 = fread(fid, 1, 'short')';
header.scalar_name               = fread(fid, [20,10], '*char')';
header.n_properties              = fread(fid, 1, 'short')';
header.property_name             = fread(fid, [20,10], '*char')';
header.reserved                  = fread(fid, 508, '*char');
header.voxel_order               = fread(fid, 4, '*char')';
header.pad2                      = fread(fid, 4, '*char')';
header.image_orientation_patient = fread(fid, 6, 'float')';
header.pad1                      = fread(fid, 2, '*char')';
header.invert_x                  = fread(fid, 1, 'uchar');
header.invert_y                  = fread(fid, 1, 'uchar');
header.invert_z                  = fread(fid, 1, 'uchar');
header.swap_xy                   = fread(fid, 1, 'uchar');
header.swap_yz                   = fread(fid, 1, 'uchar');
header.swap_zx                   = fread(fid, 1, 'uchar');
header.n_count                   = fread(fid, 1, 'int')';
header.version                   = fread(fid, 1, 'int')';
header.hdr_size                  = fread(fid, 1, 'int')';
