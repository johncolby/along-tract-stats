function [track_means,starting_pts_out,nPts] = trk_compile_data(subsDir,subIDs,tract_info,outDir,starting_pts_in, saveTrk, saveASCII)
%TRK_COMPILE_DATA - Compiles along-tract data for subjects/hemispheres/tracts
%
% Syntax: [track_means,starting_pts_out,nPts] = trk_compile_data(subsDir,subIDs,tract_info,outDir,starting_pts_in,saveTrk,saveASCII)
%
% Inputs:
%    subsDir    - Path to subject directory [char]
%    subIDs     - List of subject ID folders in subsDir [nSubs x 1 cell or num array]
%    tract_info - Dataset with tract names, default origins, and viewpoints
%    outDir     - Path to output directory. (Default: PWD)
%    starting_pts_in - Dataset with tract origins. Useful for reusing past
%        origins that were determined interactively
%    saveTrk    - Save the mean tract geometries to .trk files? [logical]
%        (Default: 0)
%    saveASCII  - Save raw streamlines to ASCII files? [logical] (Default: 0)
%
% Outputs:
%    track_means      - Structure array with mean tract geometries and along-tract
%    scalar mean estimates.
%    starting_pts_out - Dataset with tract origins used for flipping. Possibly a
%        combination of defaults stored in the tract_info dataset, and
%        interactively chosen ones.
%    nPts             - The number of interpolation points used for each tract
%                       and subject [nTracts x nSubjects]
%
% Output files:
%    starting_pts_out.txt - Same as above
%    Tracking_QC_*.pdf - Quality control images for each tract. Allows you to
%        quickly see if the flipping/interpolating appears successful.
%    trk_data.txt - Data table with scalar means and standard deviations for
%        each subject/tract/hemisphere (for stats in R).
%    trk_props_long.txt - Data table with number of streamlines for each
%        subject/tract/hemisphere (for stats in R).
%    <subID>_<trkName>.txt ? (Optionally) Raw streamlines from the original
%        tract group in a plain text ASCII format for easy plotting in R.
%        http://github.com/johncolby/along-tract-stats/wiki/single-subject
%    <subID>_<trkName>_mean.trk ? (Optionally) Mean tract geometry with attached
%        along-tract cross-sectional mean scalar metric (e.g. FA) for display in
%        TrackVis.
%        http://www.colbyimaging.com/wiki/neuroimaging/along-tract-stats#trackvis_visualization
%        
%
% Example:
%    exDir           = '/path/to/along-tract-stats/example';
%    subIDs          = {'subject1'};
%    tract_info      = dataset('file', fullfile(exDir, 'tract_info.txt'));
%    starting_pts_in = dataset('file', fullfile(exDir, 'starting_pts_out.txt'));
%    [track_means,starting_pts_out,nPts] = trk_compile_data(exDir,subIDs,tract_info,[],starting_pts_in,1,1);
%
% Other m-files required: read_avw, trk_read, trk_interp, trk_flip, trk_restruc,
% trk_add_sc, trk_mean_sc, trk_plot, trk_write_ascii, dataset, export
% Subfunctions: none
% MAT-files required: none
%
% See also: DATASET

% Author: John Colby (johncolby@ucla.edu)
% UCLA Developmental Cognitive Neuroimaging Group (Sowell Lab)
% Sept 2010

%% Check and format input arguments  
if nargin < 7 || isempty(saveASCII), saveASCII = 0; end
if nargin < 6 || isempty(saveTrk), saveTrk = 0; end
if nargin < 5 || isempty(starting_pts_in)
    starting_pts_in = dataset();
elseif ~iscell(starting_pts_in.Subject) % Reformat to cells if subIDs are all numeric
    starting_pts_in.Subject = cellfun(@num2str, num2cell(starting_pts_in.Subject), 'UniformOutput', 0);
end
if nargin < 4 || isempty(outDir), outDir = pwd; end
if nargin < 3, error('Need more input arguments'), end
if isnumeric(subIDs), subIDs = cellstr(num2str(subIDs)); end

tract_info.startPt = cell2mat(cellfun(@eval, tract_info.startPt, 'UniformOutput',false));
tract_info.view    = cell2mat(cellfun(@eval, tract_info.view, 'UniformOutput',false));

%% Setup output files and variables
% Save one file for overall track properties
fid1 = fopen(fullfile(outDir, 'trk_props_long.txt'), 'wt');
fprintf(fid1, 'ID\tHemisphere\tTract\tStreamlines');

% Save another file with the mean FA +/- SD at each point along the track
fid2 = fopen(fullfile(outDir, 'trk_data.txt'), 'wt');
fprintf(fid2, 'ID\tPoint\tHemisphere\tTract\tFA\tSD');

% Initialize variables
track_means      = struct([]);
starting_pts_out = dataset();
nPts             = zeros(length(tract_info),length(subIDs));

%% Main loop to extract along-tract properties
% Loop over tracks
for iTrk=1:length(tract_info)
    ipage   = 1;
    newpage = 0;
    fh      = figure; hold on
    
    % Loop over subjects
    for i=1:length(subIDs)
        pt_start = [];
        if (i-1)/25 == ipage, newpage = 1; end
        
        try
            % Load scalar volume
            % Note: Modify path according to your directory setup
            subStr  = subIDs{i};
            trkName = sprintf('%s_%s', tract_info.Tract{iTrk}, tract_info.Hemisphere{iTrk});
            volPath = fullfile(subsDir, subStr, 'dti_fa.nii.gz');
            volume  = read_avw(volPath);
            
            % Load tract group
            % Note: Modify path according to your directory setup
            trkPath = fullfile(subsDir, subStr, strcat(trkName, '.trk'));
            [header tracks] = trk_read(trkPath);
            
            % Determine # of interpolation points
            if isnan(tract_info.nPts(iTrk))
                nPts(iTrk,i) = round(mean(trk_length(tracks))/header.voxel_size(1));
            else
                nPts(iTrk,i) = tract_info.nPts(iTrk);
            end
            
            % Interpolate streamlines
            tracks_interp = trk_interp(tracks, nPts(iTrk,i), [], 1);
            nPts(iTrk,i)  = size(tracks_interp, 1);
            
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
            for iPt=1:nPts(iTrk,i)
                fprintf(fid2, '\n%s\t%d\t%s\t%s\t%0.4f\t%0.4f', subStr, iPt, tract_info.Hemisphere{iTrk}, tract_info.Tract{iTrk}, scalar_mean(iPt), scalar_sd(iPt));
            end
            
            % Determine the mean streamline geometry for display in QC figures
            track_mean = mean(tracks_interp, 3);
            track_mean_sc = [track_mean scalar_mean];
            track_means = [track_means struct('Subject', subStr, 'Tract', tract_info.Tract{iTrk}, 'Hemisphere', tract_info.Hemisphere{iTrk}, 'track_mean_sc_str', trk_restruc(track_mean_sc))];
            
            % Open a new QC figure if needed
            if newpage
                set(gcf, 'PaperSize', [10.5 8])
                set(gcf, 'PaperPosition', [0 0 10.5 8])
                print(gcf, '-dpdf', fullfile(outDir, sprintf('Tracking_QC_%s_%d.pdf', trkName, ipage)), '-r300')
                close(fh)
                fh = figure; hold on
                ipage = ipage+1;
                newpage = 0;
            end
            
            % Draw QC figure
            figure(fh)
            subplot(5,5,i-25*(ipage-1))
            trk_plot(header, trk_restruc(track_mean_sc), volume, [])
            view(tract_info.view(iTrk,:))
            title(subStr)
            axis off
                        
            % Save the mean tract geometry if desired
            if saveTrk
                %                         x y z   sc1
                track_mean_sc(1:2,:,2) = [0 0 0   0;  %min
                                          0 0 0.1 1]; %max
                track_mean_sc_str = trk_restruc(track_mean_sc);
                
                header_mean_sc = header_sc;
                header_mean_sc.n_count = 2;
                trk_write(header_mean_sc, track_mean_sc_str, fullfile(outDir, sprintf('%s_%s_mean.trk', subStr, trkName)))
            end
            
            % Save raw streamlines to ASCII if desired
            if saveASCII
                tracks_sc_mat = trk_restruc(tracks_sc);
                trk_write_ascii(tracks_sc_mat, fullfile(outDir, sprintf('%s_%s.txt', subStr, trkName)))
            end
            
        catch me % No streamlines
            fprintf(fid1, '\n%s\t%s\t%s\t0', subStr, tract_info.Hemisphere{iTrk}, tract_info.Tract{iTrk});
            fprintf('Failed to process subject %s %s\n', subStr, trkName)
            warning(me.message)
        end
    end
    
    % Save QC figure
    set(gcf, 'PaperSize', [10.5 8])
    set(gcf, 'PaperPosition', [0 0 10.5 8])
    print(gcf, '-dpdf', fullfile(outDir, sprintf('Tracking_QC_%s_%d.pdf', trkName, ipage)), '-r300')
    close(fh)
end

%% Clean up
export(starting_pts_out, 'file', fullfile(outDir, 'starting_pts_out.txt'))

fclose all;
