function trk_stats_overlay(efxFile,subDir,outDir,starting_pts_in)
%TRK_STATS_OVERLAY - Overlay stats on an example subject's mean tract geometry
%This function overlays the statistical results from R (effect sizes and
%p-values) onto an example subject's mean tract geometry. The output is saved as
%a .trk file for final visualization in TrackVis.
%
% Syntax: trk_stats_overlay(efxFile,subDir,outDir,starting_pts_in)
%
% Inputs:
%    efxFile - Path to 'effects_table.txt' from R
%    subDir  - Path to directory for the example subject
%    outDir  - Path to output directory (Default: PWD)
%    starting_pts_in - Dataset with tract origins. Useful for reusing past
%        origins that were determined interactively
%
% Output files:
%    results.trk - Output tract group with 1 streamline for each tract/hemisphere
%        in effects_table. Attached scalars include the effect sizes and p-values
%        Load into TrackVis to visualize. 
%
% Example:
%    exDir           = '/path/to/along-tract-stats/example';
%    efxFile         = fullfile(exDir, 'effects_table.txt');
%    subDir          = fullfile(exDir, 'subject1');
%    starting_pts_in = dataset('file', fullfile(exDir, 'starting_pts_out.txt'));
%    trk_stats_overlay(efxFile,subDir,[],starting_pts_in)
%
% Other m-files required: read_avw, trk_read, trk_interp, trk_flip, trk_restruc,
% trk_plot, trk_write, dataset
% Subfunctions: none
% MAT-files required: none
%
% See also: TRK_COMPILE_DATA

% Author: John Colby (johncolby@ucla.edu)
% UCLA Developmental Cognitive Neuroimaging Group (Sowell Lab)
% Sept 2010

%% Check and format input arguments
if nargin < 4 || isempty(starting_pts_in)
    starting_pts_in = dataset();
    warning('You should really use the same tract origins as when you processed the data. Is that ''starting_pts'' dataset available?')
elseif ~iscell(starting_pts_in.Subject) % Reformat to cells if subIDs are all numeric
    starting_pts_in.Subject = cellfun(@num2str, num2cell(starting_pts_in.Subject), 'UniformOutput', 0);
end
if nargin < 3 || isempty(outDir), outDir = pwd; end
if nargin < 2, error('Need more input arguments'), end

%% Load files and setup
% Load the statistical results dataset from R
ds          = dataset('file', efxFile, 'Delimiter', ' ');
Tracts      = unique(ds.Tract);
Hemispheres = unique(ds.Hemisphere);
vars        = {'Value' 'p0x2Evalue' 'p0x2Eval0x2Eadj'};
names       = char('effect', 'pval', 'pval_corr');
% Begin with a dummy tract to set scalar ranges
                         %x y z      sc1  sc2  sc3
track_mean_sc(1:2,:,1) = [0 0 0     -0.1  0    0;  %min
                          0 0 0.1    0.1  1    1]; %max
% Trim if a column of corrected p-values isn't present
if ~any(strcmp(get(ds, 'VarNames'), 'p0x2Eval0x2Eadj'))
    vars  = vars(1:2);
    names = names(1:2,:);
    track_mean_sc = track_mean_sc(:,1:5,1);
end
track_mean_sc_str = trk_restruc(track_mean_sc);

% Load the scalar volume
% Note: Modify path according to your directory setup
volPath = fullfile(subDir, 'dti_fa.nii.gz');
volume  = read_avw(volPath);

%% Generate mean tract geometry and overlay statistical results
% Loop over each tract and hemisphere
for iTrk=1:length(Tracts)
    for iHemi=1:length(Hemispheres)
        % Load up a representative tract group
        trkFile = sprintf('%s_%s.trk', Tracts{iTrk}, Hemispheres{iHemi});
        trkPath = fullfile(subDir, trkFile);
        [header tracks] = trk_read(trkPath);
        
        % Slice out this tract from the effects table
        effects = ds(strcmp(Tracts{iTrk}, ds.Tract) & strcmp(Hemispheres{iHemi}, ds.Hemisphere),:);
        
        % Determine how many interpolation point were used
        nPts = length(effects);
        
        % Reorient and interpolate the tract group
        tracks_interp = trk_interp(tracks, nPts);
        
        % Determine 'starting_pt' near tract origin. First look in
        % 'starting_pts_in' if available.
        pt_start = [];
        [~,subStr,~,~] = fileparts(subDir);
        if ~isempty(starting_pts_in)
            pt_start = double(starting_pts_in(strcmp(subStr, starting_pts_in.Subject) &...
                strcmp(Tracts{iTrk}, starting_pts_in.Tract) &...
                strcmp(Hemispheres{iHemi}, starting_pts_in.Hemisphere),4:6));
        end
        
        % Reorient streamlines according to 'pt_start'. Determine
        % interactively if needed
        [tracks_interp pt_start] = trk_flip(header, tracks_interp, pt_start, volume);
        
        % Extract the effect size and p-value scalars for this tract group
        data = double(effects(:,vars));
        
        % Generate mean streamline gemoetry and attach scalars
        track_mean = mean(tracks_interp, 3);
        track_mean_sc_str = [track_mean_sc_str trk_restruc([track_mean data])];
    end
end

%% Outputs
% Generate a new header with the right number of streamlines and scalars
header_mean_sc           = header;
header_mean_sc.n_count   = length(track_mean_sc_str);
header_mean_sc.n_scalars = size(names,1);
header_mean_sc.scalar_name(1:size(names, 1),1:size(names,2)) = names;

% Save the results for final visualization in TrackVis
savePath = fullfile(outDir, 'results.trk');
trk_write(header_mean_sc, track_mean_sc_str, savePath)
