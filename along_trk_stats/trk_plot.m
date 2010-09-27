function trk_plot(header,tracks,volume,slices)
%TRK_PLOT - 3D plot of TrackVis .trk track group
%
% Syntax: trk_plot(header,tracks)
%
% Inputs:
%    header - .trk file header 
%    tracks - .trk file body
%    volume - (optional) Scalar MRI volume to use for slice overlays
%    slice  - (optional) XYZ slice planes (in voxels) for overlays [1 x 3] 
%             (Default = [0 0 0])
%
% Outputs:
%
% Example: 
%    [header tracks] = trk_read(filePath)
%    trk_plot(header, tracks)
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: TRK_READ

% Author: John Colby (johncolby@ucla.edu)
% UCLA Developmental Cognitive Neuroimaging Group (Sowell Lab)
% Mar 2010 $Rev$ $Date$

% Input argument defaults
if nargin == 3, slices = [0 0 0]; end

if ~isstruct(tracks), error('Tracks must be in structure form. Try running TRK_RESTRUC first.'), end

hold on
for iTrk = 1:length(tracks)
   %plot3(tracks(iTrk).matrix(:,1)/header.voxel_size(1), tracks(iTrk).matrix(:,2)/header.voxel_size(2), tracks(iTrk).matrix(:,3)/header.voxel_size(3))
    plot3(tracks(iTrk).matrix(:,1), tracks(iTrk).matrix(:,2), tracks(iTrk).matrix(:,3))
    plot3(tracks(iTrk).matrix(1,1), tracks(iTrk).matrix(1,2), tracks(iTrk).matrix(1,3), 'r.')
end

% Plot slice overlays
if nargin>2 && ~isempty(volume)
    slices(2) = header.dim(2) - slices(2); % Temporary tweak until I can figure
    slices    = slices.*header.voxel_size; % out the orientation issues
    [x, y, z] = meshgrid(header.voxel_size(1)*(0:(header.dim(1)-1)),...
                         header.voxel_size(2)*(0:(header.dim(2)-1)),...
                         header.voxel_size(3)*(0:(header.dim(3)-1)));
    slice(x,y,z,permute(volume, [2 1 3]), slices(1), slices(2), slices(3));
    shading flat
    colormap(gray)
end

xlabel('x'), ylabel('y'), zlabel('z', 'Rotation', 0)
box on
axis image
axis ij
