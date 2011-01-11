function [track_means,scalar_means,starting_pts_out] = trk_compile_data(subsDir,subIDs,tract_info,nPts,outDir,starting_pts_in)
%TRK_COMPILE_DATA - Compiles along-tract data for subjects/hemispheres/tracts
%This function performs the along-tract
%
% Syntax: [track_means,scalar_means,starting_pts_out] = trk_compile_data(subsDir,subIDs,tract_info,nPts,outDir,starting_pts_in)
%
% Inputs:
%    subsDir    - Path to subject directory [char]
%    subIDs     - List of subject ID folders in subsDir [n x 1 cell or num array]
%    tract_info - Dataset with tract names, default origins, and viewpoints
%    nPts       - Number of points for tract interpolation (Default: 30)
%    outDir     - Path to output directory. (Default: PWD)
%    starting_pts_in - Dataset with tract origins. Useful for reusing past
%        origins that were determined interactively
%
% Outputs:
%    track_means      - 4D matrix of mean tract geometries [nPts x 3 x nTracts x 
%        nSubjects]
%    scalar_means     - 3D matrix of tract scalars [nPts x nTracts x nSubjects]
%    starting_pts_out - Tract origins used for flipping. Possibly a combination
%        of defaults stored in the tract_info dataset, and interactively chosen
%        ones.
%
% Output files:
%    starting_pts_out.txt - Same as above
%    Tracking_QC_*.pdf - Quality control images for each tract. Allows you to
%        quickly see if the flipping/interpolating appears successful.
%    trk_data.txt - 
%    trk_props_long.txt - 
%
% Example: 
%    subsDir = '/ifs/edevel/TRIO/DATA_ANALYSES/JC_MOD/ANALYSIS/prelim/SUBJECTS';
%    subIDs  = num2cell([20037; 20103]);
%    tract_info      = dataset('file', '/Users/jcolby/Documents/LONI/stats_along_tracts/tract_info.txt');
%    starting_pts_in = dataset('file', '/ifs/edevel/TRIO/DATA_ANALYSES/JC_MOD/starting_pts.txt');
%    [track_means,scalar_means,starting_pts_out] = trk_compile_data(subsDir,subIDs,tract_info,[],[],starting_pts_in);
%
% Other m-files required: read_avw, trk_read, trk_interp, trk_flip, trk_restruc,
% trk_add_sc, trk_mean_sc, trk_plot, dataset, export
% Subfunctions: none
% MAT-files required: none
%
% See also: DATASET

% Author: John Colby (johncolby@ucla.edu)
% UCLA Developmental Cognitive Neuroimaging Group (Sowell Lab)
% Sept 2010

%% Check and format input arguments
if nargin < 6 || isempty(starting_pts_in)
    starting_pts_in = dataset();
elseif ~iscell(starting_pts_in.Subject) % Reformat to cells if subIDs are all numeric
    starting_pts_in.Subject = cellfun(@num2str, num2cell(starting_pts_in.Subject), 'UniformOutput', 0);
end
if nargin < 5 || isempty(outDir), outDir = pwd; end
if nargin < 4 || isempty(nPts), nPts = 30; end
if nargin < 3, error('Need more input arguments'), end
if isnumeric(subIDs), subIDs = num2cell(subIDs); end

tract_info.startPt = cell2mat(cellfun(@eval, tract_info.startPt, 'UniformOutput',false));
tract_info.view    = cell2mat(cellfun(@eval, tract_info.view, 'UniformOutput',false));

%% Setup output files and variables
% Save one file for overall track properties
fid1 = fopen(fullfile(outDir, 'trk_props_long.txt'), 'wt');
fprintf(fid1, 'ID\tHemisphere\tTract\tStreamlines');

% Save another file with the mean FA +/- SD at each point along the track
fid2 = fopen(fullfile(outDir, 'trk_data.txt'), 'wt');
fprintf(fid2, 'ID\tPoint\tHemisphere\tTract\tFA\tSD');

% Initialize QC figures
for iFig=1:length(tract_info), fh(iFig) = figure; hold on; end

% Initialize storage variables
track_means      = zeros(nPts, 3, length(tract_info), length(subIDs));
scalar_means     = zeros(nPts, length(tract_info), length(subIDs));
starting_pts_out = dataset();

%% Main loop to extract along-tract properties
% Loop over subjects
for i=1:length(subIDs)
    subStr  = num2str(subIDs{i});
    
    % Load scalar volume
    % Note: Modify path according to your directory setup
    volPath = fullfile(subsDir, subStr, 'DTI/diffusion_toolkit/dti_fa.nii.gz');
    volume  = read_avw(volPath);
    
    % Loop over tracks
    for iTrk=1:length(tract_info)
        pt_start = [];
        try
            % Load tract group
            % Note: Modify path according to your directory setup
            trkPath = fullfile(subsDir, subStr, sprintf('%s.trk', tract_info.Name{iTrk}));
            [header tracks] = trk_read(trkPath);
            
            % Interpolate streamlines
            tracks_interp = trk_interp(tracks, nPts);
            
            % Determine 'starting_pt' near tract origin. First look in
            % 'starting_pts_in' if available, and then the 'tract_info' defaults
            if ~isempty(starting_pts_in)
                pt_start = double(starting_pts_in(strcmp(subStr, starting_pts_in.Subject) &...
                    strcmp(tract_info.Tract(iTrk), starting_pts_in.Tract) &...
                    strcmp(tract_info.Hemisphere(iTrk), starting_pts_in.Hemisphere),4:6));
            end
            if isempty(pt_start)
                pt_start = tract_info.startPt(iTrk,:);
            end
            
            % Reorient streamlines according to 'pt_start'. Determine
            % interactively if needed
            [tracks_interp pt_start] = trk_flip(header, tracks_interp, pt_start, volume);
            starting_pts_out = [starting_pts_out; dataset({subStr}, tract_info.Hemisphere(iTrk), tract_info.Tract(iTrk), pt_start,...
                'VarNames',{'Subject','Hemisphere','Tract','Point'})];
            tracks_interp_str = trk_restruc(tracks_interp);
            
            % Extract scalar values from 'volume'
            [header_sc tracks_sc] = trk_add_sc(header,tracks_interp_str,volume,'FA');
            
            % Determine the mean scalar at each cross section along the tract group
            [scalar_means(:,iTrk,i) scalar_sd(:,iTrk,i)] = trk_mean_sc(header_sc,tracks_sc);
            
            % Determine the mean streamline geometry for display in QC figures
            track_means(:,:,iTrk,i) = mean(tracks_interp, 3);
            
            % Write outputs
            fprintf(fid1, '\n%s\t%s\t%s\t%d', subStr, tract_info.Hemisphere{iTrk}, tract_info.Tract{iTrk}, header.n_count);
            for iPt=1:nPts
                fprintf(fid2, '\n%s\t%d\t%s\t%s\t%0.4f\t%0.4f', subStr, iPt, tract_info.Hemisphere{iTrk}, tract_info.Tract{iTrk}, scalar_means(iPt,iTrk,i), scalar_sd(iPt,iTrk,i));
            end
        catch me % No streamlines
            fprintf(fid1, '\n%s\t%s\t%s\t0', subStr, tract_info.Hemisphere{iTrk}, tract_info.Tract{iTrk});
            fprintf('Failed to read subject %s %s\n', subStr, tract_info.Name{iTrk})
            warning(me.message)
        end
        
        % Draw QC figures
        figure(fh(iTrk))
        subplot(5,5,i)
        trk_plot(header, trk_restruc(track_means(:,:,iTrk,i)), volume, [48 48 4])
        view(tract_info.view(iTrk,:))
        title(subStr)
        axis off
    end
end

%% Save files and clean up
for iFig=1:length(tract_info)
    figure(fh(iFig))
    set(gcf, 'PaperSize', [10.5 8])
    set(gcf, 'PaperPosition', [0 0 10.5 8])
    print(gcf, '-dpdf', fullfile(outDir, sprintf('Tracking_QC_%s.pdf', tract_info.Name{iFig})), '-r300')
end

export(starting_pts_out, 'file', fullfile(outDir, 'starting_pts_out.txt'))

close all
fclose all;
