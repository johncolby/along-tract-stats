function [track_means,starting_pts_out] = trk_compile_data(subsDir,subIDs,tract_info,outDir,starting_pts_in)
%TRK_COMPILE_DATA - Compiles along-tract data for subjects/hemispheres/tracts
%
% Syntax: [track_means,starting_pts_out] = trk_compile_data(subsDir,subIDs,tract_info,outDir,starting_pts_in)
%
% Inputs:
%    subsDir    - Path to subject directory [char]
%    subIDs     - List of subject ID folders in subsDir [nSubs x 1 cell or num array]
%    tract_info - Dataset with tract names, default origins, and viewpoints
%    outDir     - Path to output directory. (Default: PWD)
%    starting_pts_in - Dataset with tract origins. Useful for reusing past
%        origins that were determined interactively
%
% Outputs:
%    track_means      - Structure array with mean tract geometries and along-tract
%    scalar mean and standard deviation estimates.
%    starting_pts_out - Dataset with tract origins used for flipping. Possibly a
%        combination of defaults stored in the tract_info dataset, and
%        interactively chosen ones.
%
% Output files:
%    starting_pts_out.txt - Same as above
%    Tracking_QC_*.pdf - Quality control images for each tract. Allows you to
%        quickly see if the flipping/interpolating appears successful.
%    trk_data.txt - Data table with scalar means and standard deviations for
%        each subject/tract/hemisphere (for stats in R).
%    trk_props_long.txt - Data table with number of streamlines for each
%        subject/tract/hemisphere (for stats in R).
%
% Example:
%    subsDir         = exDir;
%    subIDs          = {'subject1'};
%    tract_info      = dataset('file', fullfile(exDir, 'tract_info.txt'));
%    starting_pts_in = dataset('file', fullfile(exDir, 'starting_pts_out.txt'));
%    [track_means,scalar_means,starting_pts_out] = trk_compile_data(subsDir,subIDs,tract_info,[],starting_pts_in);
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
if nargin < 5 || isempty(starting_pts_in)
    starting_pts_in = dataset();
elseif ~iscell(starting_pts_in.Subject) % Reformat to cells if subIDs are all numeric
    starting_pts_in.Subject = cellfun(@num2str, num2cell(starting_pts_in.Subject), 'UniformOutput', 0);
end
if nargin < 4 || isempty(outDir), outDir = pwd; end
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

% Initialize variables
track_means = struct([]);
starting_pts_out = dataset();

%% Main loop to extract along-tract properties
% Loop over tracks
for iTrk=1:length(tract_info)
    nPts_set = 0;
    
    % Loop over subjects
    for i=1:length(subIDs)
        pt_start = [];
        
        try
            % Load scalar volume
            % Note: Modify path according to your directory setup
            subStr  = num2str(subIDs{i});
            volPath = fullfile(subsDir, subStr, 'dti_fa.nii.gz');
            volume  = read_avw(volPath);
            
            % Load tract group
            % Note: Modify path according to your directory setup
            trkFile = sprintf('%s_%s.trk', tract_info.Tract{iTrk}, tract_info.Hemisphere{iTrk});
            trkPath = fullfile(subsDir, subStr, trkFile);
            [header tracks] = trk_read(trkPath);
            
            % Determine # of interpolation points
            if nPts_set == 0
                if isnan(tract_info.nPts(iTrk))
                    nPts = round(mean(trk_length(tracks))/header.voxel_size(1));
                else
                    nPts = tract_info.nPts(iTrk);
                end
                nPts_set = 1;
            end
            
            % Interpolate streamlines
            tracks_interp = trk_interp(tracks, nPts);
            
            % Determine 'pt_start' near tract origin. First look in
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
            starting_pts_out = [starting_pts_out; dataset({subStr}, tract_info.Hemisphere(iTrk), tract_info.Tract(iTrk), pt_start(1), pt_start(2), pt_start(3),...
                'VarNames',{'Subject','Hemisphere','Tract','PointX','PointY','PointZ'})];
            tracks_interp_str = trk_restruc(tracks_interp);
            
            % Extract scalar values from 'volume'
            [header_sc tracks_sc] = trk_add_sc(header,tracks_interp_str,volume,'FA');
            
            % Determine the mean scalar at each cross section along the tract group
            [scalar_mean scalar_sd] = trk_mean_sc(header_sc,tracks_sc);
            
            % Write outputs
            fprintf(fid1, '\n%s\t%s\t%s\t%d', subStr, tract_info.Hemisphere{iTrk}, tract_info.Tract{iTrk}, header.n_count);
            for iPt=1:nPts
                fprintf(fid2, '\n%s\t%d\t%s\t%s\t%0.4f\t%0.4f', subStr, iPt, tract_info.Hemisphere{iTrk}, tract_info.Tract{iTrk}, scalar_mean(iPt), scalar_sd(iPt));
            end
            
            % Determine the mean streamline geometry for display in QC figures
            track_mean = mean(tracks_interp, 3);
            track_mean_sc_str = trk_restruc([track_mean scalar_mean scalar_sd]);
            
            % Draw QC figure
            figure
            subplot(5,5,i)
            trk_plot(header, track_mean_sc_str, volume, [])
            view(tract_info.view(iTrk,:))
            title(subStr)
            axis off
            
            % Save QC figure
            set(gcf, 'PaperSize', [10.5 8])
            set(gcf, 'PaperPosition', [0 0 10.5 8])
            trkName = sprintf('%s_%s', tract_info.Tract{iTrk}, tract_info.Hemisphere{iTrk});
            print(gcf, '-dpdf', fullfile(outDir, sprintf('Tracking_QC_%s.pdf', trkName)), '-r300')
            
        catch me % No streamlines
            fprintf(fid1, '\n%s\t%s\t%s\t0', subStr, tract_info.Hemisphere{iTrk}, tract_info.Tract{iTrk});
            fprintf('Failed to process subject %s %s\n', subStr, trkFile)
            warning(me.message)
        end
    end
end

%% Clean up
export(starting_pts_out, 'file', fullfile(outDir, 'starting_pts_out.txt'))

close all
fclose all;
