function [header tracks] = readTrk(filePath)
%READTRK   Load TrackVis .trk files

% Author: John Colby (johncolby@ucla.edu)
% Date  : 4/1/10


fid = fopen(filePath, 'r');

%% Parse in header
header.id_string                 = fread(fid, 6, '*char')';
header.dim                       = fread(fid, 3, 'short')';
header.voxel_size                = fread(fid, 3, 'float')';
header.origin                    = fread(fid, 3, 'float')';
header.n_scalars                 = fread(fid, 1, 'short')';
header.scalar_name               = fread(fid, [10,20], '*char');
header.n_properties              = fread(fid, 1, 'short')';
header.property_name             = fread(fid, [10,20], '*char');
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

%% Parse in body
tracks(header.n_count).nPoints = 0;

for iTrk = 1:header.n_count
    tracks(iTrk).nPoints = fread(fid, 1, 'int');
    tracks(iTrk).matrix  = fread(fid, [3+header.n_scalars, tracks(iTrk).nPoints], 'float')';
    if header.n_properties
        tracks(iTrk).props = fread(fid, header.n_properties, 'float');
    end
end

figure
hold on
for iTrk = 1:header.n_count
    %plot3(tracks(iTrk).matrix(:,1), tracks(iTrk).matrix(:,2), tracks(iTrk).matrix(:,3))
    plot3(tracks(iTrk).matrix(:,1)/header.voxel_size(1), tracks(iTrk).matrix(:,2)/header.voxel_size(2), tracks(iTrk).matrix(:,3)/header.voxel_size(3))

end

fclose(fid);
