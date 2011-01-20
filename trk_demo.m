%% TrackVis/MATLAB Integration Demo
% This demo takes you through the basic process of working with TrackVis track
% group .trk files in MATLAB. The main rationale for these tools is to be able
% to examine a scalar metric (e.g. FA) parameterized _along_ a track, instead
% of the typical method of collpasing across the whole track.
%
% |Author: John Colby (johncolby@ucla.edu)|
%
% |UCLA Developmental Cognitive Neuroimaging Group (Sowell Lab)|

%% Import tracks
% Set paths to .trk file and scalar MRI volume (e.g. FA map).
subDir  = fullfile(exDir, 'subject1');
trkPath = fullfile(subDir, 'CST_L.trk');
volPath = fullfile(subDir, 'dti_fa.nii.gz');

%%
% Read in the scalar MRI volume with FSL tools. (See |$FSLDIR/etc/matlab|)
volume = read_avw(volPath);

%%
% Read in the .trk file. You will get a |header| structure with fields from the .trk
% header, and a |tracks| structure array with an entry for each streamline in
% the track group.
[header tracks] = trk_read(trkPath)

%%
% Each structure in the |tracks| array contains fields for the number of points
% in that particular streamline, and the raw matrix of those points in mm 
% coordinates [nPoints x 3]. Coordinate system begins with RPI = 0,0,0.
tracks(1)
tracks(1).matrix

%% Whole-track statistics
% Determine traditional statistics collapsed across the whole track.
[meanInt stdInt nVox] = trk_stats(header, tracks, volume, 'nearest')

%%
% Determine the distribution of streamline lengths.
lengths = trk_length(tracks);
mean(lengths)
std(lengths)

%% Fit curves
% Use cubic spline interpolation to make all tracks have the same number of
% vertices. The output is in matrix mode.
tracks_interp = trk_interp(tracks, 100);

%%
% Reformat the tracks back into a structure form.
tracks_interp_str = trk_restruc(tracks_interp);

%%
% Plot the track group. The last input argument indicates the slice planes for 
% volume overlays. Red markers indicate starting points. Note the random
% ordering of tracks, with some "starting" in the cortex and some in the brainstem.
figure, trk_plot(header, tracks_interp_str, volume, [95 78 4])
view([30 30])

%%
% Reorder the tracks so that they all "start" in the same spot. The last input
% argument here is a coordinate near the desired starting area. If left blank,
% the origin will be determined interactively.
tracks_interp     = trk_flip(header, tracks_interp, [97 110 4]);
tracks_interp_str = trk_restruc(tracks_interp);

%%
% Plot the results again to see the difference.
figure, trk_plot(header, tracks_interp_str, volume, [95 78 4])
view([30 30])

%% Extract scalars
% For each point along the streamlines in |tracks_interp_str|, extract the corresponding voxel
% intensity from |volume|.
[header_sc tracks_sc] = trk_add_sc(header, tracks_interp_str, volume, 'FA');

%%
% Obtain the mean scalar at different points along the track
[scalar_mean scalar_sd] = trk_mean_sc(header_sc, tracks_sc);

%% Plot results
% Plot the mean scalar value (+/- SD) along the track.
figure, hold on
plot(scalar_mean, 'k')
plot(scalar_mean+scalar_sd, 'k--')
plot(scalar_mean-scalar_sd, 'k--')

% Make the plot prettier
hold off, box off
xlim([0 100]), ylim([0 1])
title('\bf{Mean FA along track}')
xlabel('Distance along track (%)')
ylabel('Fractional Anisotropy (FA)')

%% Export tracks
% The mean streamline geometry can also be calculated.
track_mean = mean(tracks_interp, 3);

%%
% This can be packaged up with the mean scalar values. Multiple scalars can be
% added, and can even include statistical significance values. A dummy
% streamline should be added with the desired display ranges for the scalars, if
% these exceed the ranges in the data.
track_mean_sc     = [track_mean scalar_mean];
%                         x y z   sc1
track_mean_sc(1:2,:,2) = [0 0 0   0;  %min
                          0 0 0.1 1]; %max
track_mean_sc_str = trk_restruc(track_mean_sc);

%%
% Don't forget to update the header because now we only have 1 streamline
% instead of many (plus the 1 dummy streamline).
header_mean_sc         = header_sc;
header_mean_sc.n_count = 2;

%%
% Save the result back out to a .trk file for visualization in TrackVis.
savePath = fullfile(exDir, 'CST_L_mean.trk');
trk_write(header_mean_sc, track_mean_sc_str, savePath)
%%
% <<mean_track.png>>


%% Help
% To get more help for any of these functions, simply type |help function_name|
% at the command prompt.
help trk_stats