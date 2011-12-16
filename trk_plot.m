function trk_plot(header,tracks,volume,slices,plottype,scalar,debug)
%TRK_PLOT - 3D plot of TrackVis .trk track group
%
% Syntax: trk_plot(header,tracks,volume,slices,plottype,scalar,debug)
%
% Inputs:
%    header  - .trk file header 
%    tracks  - .trk file body
%    volume  - (optional) Scalar MRI volume to use for slice overlays
%    slice   - (optional) XYZ slice planes (in voxels) for overlays [1 x 3] 
%              (Default: header.dim/2) Note: Use MATLAB coordinates.
%    plottype - (optional) Specify alternative plotting style. (Default is
%               empty, which simply highlights streamline origins in red)
%               [string]
%      rainbow   - Color encodes assumed correspondence, so like colors will be
%                  collapsed together.
%      scalar    - Color encode one of the tract scalars
%    scalar   - (optional) If plottype is 'scalar', specifies which scalar to
%               use (i.e. which row in header.scalar_name) [num]
%    debug    - (optional) Highlight individual vertices
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
% See also: TRK_READ, CLINE
%           http://github.com/johncolby/along-tract-stats/wiki/correspondence

% Author: John Colby (johncolby@ucla.edu)
% UCLA Developmental Cognitive Neuroimaging Group (Sowell Lab)
% Mar 2010

% Input argument defaults
if nargin < 7, debug = []; end
if nargin < 6, scalar = []; end
if nargin < 5, plottype = []; end
if nargin < 4  || isempty(slices), slices = header.dim/2; end

if ~isstruct(tracks), error('Tracks must be in structure form. Try running TRK_RESTRUC first.'), end

% Plot streamlines
hold on
maxpts = max(arrayfun(@(x) size(x.matrix, 1), tracks));
for iTrk = 1:length(tracks)
    matrix = tracks(iTrk).matrix;
    matrix(any(isnan(matrix(:,1:3)),2),:) = [];
    
    if strcmp(plottype, 'rainbow')
        cline(matrix(:,1), matrix(:,2), matrix(:,3), (0:(size(matrix, 1)-1))/(maxpts))
    elseif strcmp(plottype, 'scalar')
        matrix(isnan(matrix(:,3+scalar)),:) = [];
        cline(matrix(:,1), matrix(:,2), matrix(:,3), matrix(:,3+scalar))
    else
        plot3(matrix(:,1), matrix(:,2), matrix(:,3))
        if debug, plot3(matrix(:,1), matrix(:,2), matrix(:,3), 'r.'), end
        plot3(matrix(1,1), matrix(1,2), matrix(1,3), 'r.')
    end
end

% Plot slice overlays
if nargin>2 && ~isempty(volume)
    slices    = (slices - 0.5).*header.voxel_size;
    [x, y, z] = meshgrid(header.voxel_size(1)*(0.5:header.dim(1)),...
                         header.voxel_size(2)*(0.5:header.dim(2)),...
                         header.voxel_size(3)*(0.5:header.dim(3)));
    h2 = slice(x,y,z,permute(volume, [2 1 3]), slices(1), slices(2), slices(3), 'nearest');
    shading flat
    if any(strcmp(plottype, {'rainbow' 'scalar'}))
        colormap([jet(100);gray(100)])
        slice_cdata = get(h2, 'CData');
        slice_cdata = cellfun(@(x) x+1, slice_cdata, 'UniformOutput', false);
        for iSlice=1:3
            set(h2(iSlice), 'CData', slice_cdata{iSlice});
        end
        caxis([0 2])
    else
        colormap(gray)
    end
end

xlabel('x'), ylabel('y'), zlabel('z', 'Rotation', 0)
box on
axis image
axis ij
