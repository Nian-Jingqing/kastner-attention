%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%LIP_Exo_target_cued: Target locked (-1000:500) LFP at RFin location
%lipmat: LIP channels
% LIP_Exo_response: array of responses of trials (1 hit 0 miss) 
%frange: freqeuncey range of interest
%srate: sampling rate
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% For one session
PhaseBehav_tmp1=[];
for j1=1:16
    tmp1(1,:,:)=LIP_Exo_target_cued{lipmat(j1)};
    [fp1,famp1]=get_Hilbert(tmp1,srate,frange);
    tmp2=[rad2deg(wrapTo2Pi(squeeze(fp1(:,1001,:)))) LIP_Exo_response{lipmat(j1)}];
    PhaseBehav_tmp1=[PhaseBehav_tmp1;tmp2];tmp2=[];
    tmp1=[];fp1=[];famp1=[];
end
LIP_PhaseBehav{ik}=PhaseBehav_tmp1;PhaseBehav_tmp1=[];
%%%%%%%%%%%%%%%%%%%%%
ii1=0:5:180;
ii2=wrapTo360(ii1+180);
ii3=[ii1' ii2'];
ii4=180:5:360;
ii5=wrapTo360(ii4+180);
ii6=[ii4' ii5'];
%PhaseDetection - LIP
for jj=1:37
    resp_tmp1=LIP_PhaseBehav(:,2);
    resp_tmp2=resp_tmp1(LIP_PhaseBehav(:,1) >= ii3(jj,1) & LIP_PhaseBehav(:,1) <= ii3(jj,2));
    LIP_PhaseDetection1(jj)=sum(resp_tmp2(resp_tmp2==1))/length(resp_tmp2);
    resp_tmp1=[];resp_tmp2=[];
end
for jj=1:37
    resp_tmp1=LIP_PhaseBehav(:,2);
    resp_tmp2=resp_tmp1(LIP_PhaseBehav(:,1) >= ii6(jj,1) | LIP_PhaseBehav(:,1) <= ii6(jj,2));
    LIP_PhaseDetection2(jj)=sum(resp_tmp2(resp_tmp2==1))/length(resp_tmp2);
    resp_tmp1=[];resp_tmp2=[];
end
LIP_PhaseDetection_tmp1=([LIP_PhaseDetection1 LIP_PhaseDetection2]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
subplot(221)
f = fit((1:74)',LIP_PhaseDetection_tmp2', 'sin2')
plot(f,(1:74)',LIP_PhaseDetection_tmp2')
ax=gca;
ax.XTick = [1 18 36 54 72];ax.XTickLabel={'0','90','180','270','360'};xlim([1 74]);
ylabel('Hit Rate')
xlabel('ThetaPhase')
set(gca, 'FontSize', 18);
set(gcf, 'DefaultAxesFontName', 'Arial');
set(gcf, 'DefaultTextFontName', 'Arial');
set(gca, 'box','off')
set(gcf,'color','w');h=gcf;
set(h,'PaperOrientation','landscape');
set(h,'PaperUnits','normalized');
set(h,'PaperPosition', [0 0 1 1]);
print('-depsc','PhaseDetection-LIP');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
