function trk_stats_overlay(exptDir,subDir,tract_info,corrected)
%TRK_STATS_OVERLAY - Overlay stats on an example subject's mean tract geometry
%This function overlays the statistical results from R (effect sizes and
%p-values) onto an example subject's mean tract geometry. The output is saved as
%a .trk file for final visualization in TrackVis.
%
% Syntax: trk_stats_overlay(exptDir,subDir,tract_info)
%
% Inputs:
%    exptDir    - Path to directory with 'effects_table.txt' file from R
%    subDir     - Path to directory for the example subject
%    tract_info - Dataset with tract names, default origins, and viewpoints
%
% Output files:
%    results.trk - Saved in 'exptDir'. Load into TrackVis to visualize. 
%
% Example: 
%    exptDir    = '/Users/jcolby/Documents/LONI/stats_along_tracts';
%    subDir     = '/ifs/edevel/TRIO/DATA_ANALYSES/JC_MOD/SUBJECTS/20732';
%    tract_info = dataset('file', '/Users/jcolby/Documents/LONI/stats_along_tracts/tract_info.txt');
%    trk_stats_overlay(exptDir,subDir,tract_info)
%
% Other m-files required: read_avw, trk_read, trk_interp, trk_flip, trk_restruc,
% trk_plot, trk_write, dataset
% Subfunctions: none
% MAT-files required: none
%
% See also: TRK_COMPILE_DATA

% Author: John Colby (johncolby@ucla.edu)
% UCLA Developmental Cognitive Neuroimaging Group (Sowell Lab)
% Sept 2010 $Rev$ $Date$

%% Load files and setup
% Load the statistical results dataset from R
ds = dataset('file', fullfile(exptDir, 'effects_table.txt'), 'Delimiter', ' ');
Tracts      = unique(ds.Tract);
Hemispheres = unique(ds.Hemisphere);

% Load the scalar volume
% Note: Modify path according to your directory setup
volPath = fullfile(subDir, 'DTI/diffusion_toolkit/dti_fa.nii.gz');
volume  = read_avw(volPath);

% Begin with a dummy tract to set scalar ranges
track_mean_sc = zeros(30,5);
                         %x y z      sc1  sc2
track_mean_sc(1:2,:,1) = [0 0 0     -0.1  0;  %min
                          0 0 0.1    0.1  1]; %max

%% Generate mean tract geometry and overlay statistical results
% Loop over each tract and hemisphere
for iTrk=1:length(Tracts)
    for iHemi=1:length(Hemispheres)
        % Load up a representative tract group
        trkFile = sprintf('%s.trk', tract_info.Name{strcmp(Hemispheres{iHemi}, tract_info.Hemisphere) & strcmp(Tracts{iTrk}, tract_info.Tract)});
        trkPath = fullfile(subDir, trkFile);
        [header tracks] = trk_read(trkPath);
        
        % Reorient and interpolate the tract group
        tracks_interp     = trk_interp(tracks, 30);
        tracks_interp     = trk_flip(header, tracks_interp, [], volume);
        tracks_interp_str = trk_restruc(tracks_interp);
        
        trk_plot(header, tracks_interp_str, volume, [48 48 1])
        
        % Extract the effect size and p-value scalars for this tract group
        if corrected==1
            data = double(ds(strcmp(Tracts{iTrk}, ds.Tract) & strcmp(Hemispheres{iHemi}, ds.Hemisphere),{'Value' 'p0x2Eval0x2Eadj'}));
        else
            data = double(ds(strcmp(Tracts{iTrk}, ds.Tract) & strcmp(Hemispheres{iHemi}, ds.Hemisphere),{'Value' 'p0x2Evalue'}));
        end
        
        % Generate mean streamline gemoetry and attach scalars
        track_mean    = mean(tracks_interp, 3);
        track_mean_sc = cat(3, track_mean_sc, [track_mean data]);
    end
end

track_mean_sc_str = trk_restruc(track_mean_sc);

%% Outputs
% Generate a new header with the right number of streamlines and scalars
header_mean_sc           = header;
header_mean_sc.n_count   = length(track_mean_sc_str);
header_mean_sc.n_scalars = 2;
name = 'effect';
header_mean_sc.scalar_name(1,1:size(name,2)) = name;
name = 'pval';
header_mean_sc.scalar_name(2,1:size(name,2)) = name;

% Save the results for final visualization in TrackVis
savePath = fullfile(exptDir, 'results.trk');
trk_write(header_mean_sc, track_mean_sc_str, savePath)