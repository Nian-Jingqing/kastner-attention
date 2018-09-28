function [] = generateTrialPlot(R,trialNum)
    trialIdx = trialNum * 2;
    Rx = [R(trialIdx).x];
    Rtarg =  [R(trialIdx).target];
    Remg =  [R(trialIdx).emg];
    length = 1:size(Rx,2);
    figure()
    subplot(2,3,1)
    plot(length,Rtarg(1,:))
    hold on
    plot(length,Rx(1,:))
    plot(length,Rx(3,:))
    hold off
    title('X Position')
    legend('Target Pos.','Monkey Pos.','Monkey Vel.','Location','Best')
    subplot(2,3,2)
    plot(length,Rtarg(2,:))
    hold on
    plot(length,Rx(2,:))
    plot(length,Rx(4,:))
    hold off
    title('Y Position')
    legend('Target Pos.','Monkey Pos.','Monkey Vel.','Location','Best')    
    subplot(2,3,3)
    plot(length,Remg(1,:))
    hold on
    plot(length,Remg(2,:))
    hold off
    title('EMG Activation- cPec latD')
    legend('cPec','latD','Location','Best')
    subplot(2,3,4)
    plot(length,Remg(3,:))
    hold on
    plot(length,Remg(4,:))
    hold off
    title('EMG Activation- aDel mDel')
    legend('aDel','mDel','Location','Best')
    subplot(2,3,5)
    plot(length,Remg(5,:))
    hold on
    plot(length,Remg(6,:))
    plot(length,Remg(7,:))
    hold off
    title('EMG Activation- mBic lBic lTri')
    legend('mBic','lBic','lTri','Location','Best')
    subplot(2,3,6)
    hold on
    plot(length,Remg(8,:))
    plot(length,Remg(9,:))
    plot(length,Remg(10,:))
    hold off
    title('EMG Activation- brac exCR flCR')
    legend('brac','exCR','flCR','Location','Best')
end
