%%
%Projection;
%FRandSpiking;
alignPoint = 200/binsize_rescaled;
windowStart = 520/binsize_rescaled;
windowEnd = 800/binsize_rescaled;

savedir3 = '/snel/share/share/derived/kastner/LFADS_runs/pulvinar/Multi-day/multiDay_CO_AO_TD_HoldRel_JanToApr/postAnalysis/withExternalInput_20180614/jPCA_Rasters/tmp';
%%
for nt = 1 : 5 : numel( Projection )
    %for nt = 2
    f1 = figure
    sp1 = subplot(4,1,1);
    plot( Projection( nt ).projAllTimes(:,3)');
    set(gca,'XTick',[1 alignPoint windowStart windowEnd]);
    set(gca,'XTickLabels',{'-200ms', 'ArrayOn', '+320ms','+600ms'});
    ylabel('Dim3')
    title(sp1, 'jPC3')
    hold on

    sp2 = subplot(4,1,2);
    plot( Projection( nt ).projAllTimes(:,4)');
    set(gca,'XTick',[1 alignPoint windowStart windowEnd]);
    set(gca,'XTickLabels',{'-200ms', 'ArrayOn', '+320ms','+600ms'});
    ylabel('Dim4')
    title(sp2, 'jPC4')

    norm_FR = normalize(FRandSpiking(nt).FR, 'centered');
    
    sp3 = subplot(4, 1, 3);
    imagesc( norm_FR );
    set(gca,'XTick',[1 alignPoint windowStart windowEnd]);
    set(gca,'XTickLabels',{'-200ms', 'ArrayOn', '+320ms','+600ms'});
    title(sp3, 'LFADS rates');
    ylabel('Multi-units');

    sp4 = subplot(4, 1, 4);
    imagesc( FRandSpiking(nt).spiking );
    set(gca,'XTick',[1 alignPoint windowStart windowEnd]);
    set(gca,'XTickLabels',{'-200ms', 'ArrayOn', '+320ms','+600ms'});
    ylabel('Multi-units');
    title(sp4, 'Real Spiking');

    suptitle(['Trial ' int2str(nt)])
    set(f1, 'Position', [310 4 1142 962]);
    cd(savedir3)
    print(f1,['Trial ' int2str(nt)], '-dpng');

    close;
end
