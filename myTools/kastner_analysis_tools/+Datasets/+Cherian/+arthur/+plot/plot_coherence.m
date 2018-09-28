function [ ] = plot_coherence( f, cxy, channel, preds, freq_range )

figure()
for i = 1:numel( preds )
    i_pred = preds( i );
    plot( f, squeeze( cxy( i_pred, channel, : ) ) )
    hold on
end
hold off
xlim( freq_range )
xlabel( 'Frequency (Hz)' )
ylabel( 'Magnitude-Squared Coherence' )

