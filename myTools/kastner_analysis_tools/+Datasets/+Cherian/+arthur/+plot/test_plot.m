function [ ] = test_plot( output, r_succ, full_neural, new_sr, figs_dir, i_trial )
addpath( '/snel/home/lwimala/Projects/MATLAB/dependencies/subaxis')
orig_sr = r_succ.r( i_trial ).sampleRate;
orig_start_idx = r_succ.r( i_trial ).abstSidx;
orig_end_idx = r_succ.r( i_trial ).abstEidx;

t1_idx = r_succ.r( i_trial ).relt1idx;
t2_idx = r_succ.r( i_trial ).relt2idx;
t3_idx = r_succ.r( i_trial ).relt3idx;
tE_idx = r_succ.r( i_trial ).reltEidx;

%disp( t1_idx / tE_idx )
%disp( t2_idx / tE_idx )
%disp( t3_idx / tE_idx )
%disp( 500 / tE_idx )

disp( r_succ.r( i_trial ).trialNum )


start_idx = round( orig_start_idx * ( new_sr / orig_sr ) );
end_idx = round( orig_end_idx * ( new_sr / orig_sr ) );

plot_range = start_idx:end_idx;

%% plot neural

figSizeX = 10;
figSizeY = 10;
figure_prop_name = { 'PaperPositionMode','units','Position' };
figure_prop_val = { 'auto', 'inches', [ 1 1 figSizeX figSizeY ] };


set( gcf, figure_prop_name, figure_prop_val )
axes('Units', 'normalized', 'Position', [0 0 1 1])

for i = 1:numel( full_neural )
    vertSpace = 0.0;
    % plot binned spikes
    if i == 1
        data = full_neural{ i }( plot_range, : );
        bool_data = ( data > 0 );
        bool_data = ~bool_data;
        
        subaxis( numel( full_neural ), 1, i, 'SpacingVert', vertSpace)
        imagesc( bool_data' )
        colormap( bone )
    else
        data = full_neural{ i }( plot_range, : );
        if i == 2
            data( :, 12 ) = [];
        end
        norm_data = ( data - mean( data ) )./std( data );
        norm_data = norm_data* -1;
        subaxis( numel( full_neural ), 1, i, 'SpacingVert', vertSpace)
        imagesc( norm_data' )
        colormap( bone )
    end
    box('off')
    set( gca, 'Xtick', [ ] )
    set( gca, 'Ytick', [ ] )
    hold on
end
printFigure( figs_dir, [ '180425-arthur_robot035_s-trial_', num2str( i_trial ), '-neural' ] )

%{

tStartColor = [ 28 158 18 ]; % green
targColor1 = [ 228 191 46 ]; % yellow-gold
targColor2 = [ 225 105 20 ]; % orange
targColor3 = [ 171 11 11 ];  % deep red
targ_colors = { targColor1, targColor2, targColor3 };

% color for rates
rColor =  [ 51 200 255 ]./255; % bright blue

% color for smoothed spikes
sColor = [ 255 155 0 ]./255; % bright orange

quivColor = [ 204 251 255 ]./255;
% define colors for predictions
colors = { sColor, rColor };

opacity = { 0.8, 0.8 };
% plot position

%% initialize figure
h = Plot.blankFigure( 1 );
figSizeX = 10;
figSizeY = 10;
figure_prop_name = { 'PaperPositionMode','units','Position' };
figure_prop_val = { 'auto', 'inches', [ 1 1 figSizeX figSizeY ] };


set( gcf, figure_prop_name, figure_prop_val )
axes('Units', 'normalized', 'Position', [0 0 1 1])
xlim( [ -15 15 ] )
ylim( [ -15 15 ] )

%% plot targets

% initial x position for trial
x_0 = output{ 1 }.full_true( start_idx, 11 );
% intial y position for trial
y_0 = output{ 1 }.full_true( start_idx, 12 );

% plot initial position
arthur.plot.plotTarget( x_0, y_0, tStartColor );
hold on

for i = 1:size( r_succ.r( i_trial ).targPos, 1 )
    targX = r_succ.r( i_trial ).targPos( i, 1 );
    targY = r_succ.r( i_trial ).targPos( i, 2 );
    arthur.plot.plotTarget( targX, targY, targ_colors{ i } );
end

%% plot position

pos_x = output{ 1 }.full_true( plot_range, 11 );
pos_y = output{ 1 }.full_true( plot_range, 12 );
vel_x = output{ 1 }.full_true( plot_range, 13 );
vel_y = output{ 1 }.full_true( plot_range, 14 );

step = 20;

Plot.patchline( pos_x, pos_y, ...
                'LineWidth', 7, ...
                'edgeColor', [ 0 0 0 ], ...
                'edgeAlpha', 0.65, ...
                'faceAlpha', 0.0, ...
                'MarkerSize', 8, ...
                'MarkerFaceColor', 'k', ...
                'MarkerEdgeColor', 'k', ...
                'AlignVertexCenters', 'off' )
quiver( pos_x( 1:step:end ), pos_y( 1:step:end ), vel_x( 1:step:end ), vel_y( 1:step:end ), 'color', quivColor, ...
        'LineWidth', 4, ...
        'ShowArrowHead', 'on', ...
        'MaxHeadSize', 0.5, ...
        'AutoScaleFactor', 0.5 )

set( gca, 'visible', 'off' )    
set( gca, 'Xtick', [ ] )
set( gca, 'Ytick', [ ] )


%printFigure( figs_dir, [ '180425-arthur_robot035_s-trial_', num2str( i_trial ), '-pos' ] )
%% plot emg

figSizeX = 10;
figSizeY = 5;
figure_prop_name = { 'PaperPositionMode','units','Position' };
figure_prop_val = { 'auto', 'inches', [ 1 1 figSizeX figSizeY ] };


opacity = 0.7;
emg_name = { 'cPec', 'latD', 'aDel', 'mDel','mBic', 'lBic', 'lTri', 'brac', 'exCR', 'flCR' };
plot_channels = [ 1 2 3 6 7 8 9 ];
for i = 1:numel( plot_channels )
    i_emg = plot_channels( i );
    Plot.blankFigure( 1 + i );
    set( gcf, figure_prop_name, figure_prop_val )
    axes('Units', 'normalized', 'Position', [0 0 1 1])
    x = 1:numel( plot_range );
    for j_pred = 1: numel( output )
        if j_pred == 1
            y = output{ j_pred }.full_true( plot_range, i_emg );
            Plot.patchline( x, y, ...
                            'LineWidth', 4, ...
                            'edgeColor', [ 0 0 0 ], ...
                            'faceAlpha', 1.0, ...
                            'MarkerSize', 8, ...
                            'MarkerFaceColor', 'k', ...
                            'MarkerEdgeColor', 'k', ...
                            'AlignVertexCenters', 'off' )
            hold on
            y = output{ j_pred }.full_pred( plot_range, i_emg );
            Plot.patchline( x, y, ...
                            'LineWidth', 4, ...
                            'edgeAlpha', opacity, ...
                            'edgeColor', colors{ j_pred } , ...
                            'faceAlpha', 0.0, ...
                            'faceColor', colors{ j_pred }, ...
                            'AlignVertexCenters', 'off' )
            
        else
            hold on
            y = output{ j_pred }.full_pred( plot_range, i_emg );
            Plot.patchline( x, y, ...
                            'LineWidth', 4, ...
                            'edgeAlpha', opacity, ...
                            'edgeColor', colors{ j_pred } , ...
                            'faceAlpha', 0.0, ...
                            'faceColor', colors{ j_pred }, ...
                            'AlignVertexCenters', 'off' )
        end
    end
    set( gca, 'visible', 'off' )    
    set( gca, 'Xtick', [ ] )
    set( gca, 'Ytick', [ ] )
    %printFigure( figs_dir, [ '180425-arthur_robot035_s-trial_', num2str( i_trial ), '-pred-', emg_name{ i_emg } ] )
end
Plot.blankFigure( 7 );
%% plot neural

figSizeX = 10;
figSizeY = 10;
figure_prop_name = { 'PaperPositionMode','units','Position' };
figure_prop_val = { 'auto', 'inches', [ 1 1 figSizeX figSizeY ] };


set( gcf, figure_prop_name, figure_prop_val )
axes('Units', 'normalized', 'Position', [0 0 1 1])

for i = 1:numel( neural )
    vertSpace = 0.0;
    % plot binned spikes
    subaxis( numel( neural ), 1, i, 'SpacingVert', vertSpace)
    imagesc( neural{ i }( plot_range, : )' );
    box('off')
    set( gca, 'Xtick', [ ] )
    set( gca, 'Ytick', [ ] )
    hold on
end
%printFigure( figs_dir, [ '180425-arthur_robot035_s-trial_', num2str( i_trial ), '-neural' ] )

%}
