function [InputData] = getSpectroStims(Spectro, Stim_local,x_Average, RepStim, Flow, Win,MeanSubstractSpec)
Figure_switch=0; % Set this to one to see debugging figures

NbStim_local = length(Stim_local);
Dt = sum((1000.*Spectro.to{Stim_local(1)})<= Win);
Df = sum(Spectro.fo{Stim_local(1)}<= Flow);
x= nan(NbStim_local,Df*Dt);%this matrix will contain the vectors of spectrograms for all the stims for that window size
for ss = 1:NbStim_local
        dd=Stim_local(ss);
        %new spectro
        MatSpec = reshape(Spectro.spec{dd}, length(Spectro.fo{dd}), length(Spectro.to{dd}));
        FreqBelowFlow = find(Spectro.fo{dd}<=Flow);
        EndFreq = FreqBelowFlow(end);
        NFreq=length(FreqBelowFlow);
        TimeBelowWin = find((1000.*Spectro.to{dd})<= Win);
        EndTime = TimeBelowWin(end);
        NTime = length(TimeBelowWin);
        Newspectro=MatSpec(1:EndFreq,1:EndTime);
        if Figure_switch
            figure(100)
            imagesc(Spectro.to{dd}(1:EndTime).*1000,Spectro.fo{dd}(1:EndFreq), Newspectro)
            axis xy
            xlabel('Time (ms)')
            ylabel('Frequencies')
            title('Stim before log')
            colorbar()
            pause()
        end
        x(ss,:)=reshape(Newspectro, 1, NFreq*NTime);


end
Logx = 20*log10(abs(x));
MAXI = max(max(Logx));
Logx(Logx<(MAXI-80))=MAXI-80;
InputData.x_mean = mean(Logx);
InputData.x_std = std(Logx);
x_ZC = (Logx-repmat(InputData.x_mean,NbStim_local,1))./repmat(InputData.x_std,NbStim_local,1);
x_ZC(:,InputData.x_std==0)=0;%set these parameters to zero as they do not cary any information and would give Inf or NaN
if MeanSubstractSpec
    AvSpec_local = repmat(reshape(x_Average(1:Df, 1:Dt), 1, Df*Dt),NbStim_local,1);
    X_minusMean = Logx-AvSpec_local;
    %X_minusMean = Logx-x_Average;
    InputData.x_minusMean = nan(sum(RepStim),size(Logx,2));
end
InputData.x = nan(sum(RepStim),size(Logx,2));
InputData.x_ZC = nan(sum(RepStim),size(Logx,2));

nn=0;
for ss = 1:NbStim_local
    rr_tot = RepStim(ss);
    for rr=1:rr_tot
        nn=nn+1;
        InputData.x(nn,:) = Logx(ss,:);
        InputData.x_ZC(nn,:) = x_ZC(ss,:);
        if MeanSubstractSpec
            InputData.x_minusMean(nn,:) = X_minusMean(ss,:);
        end
    end
    
end
figure(101)
subplot(1,4,1)
imagesc(Logx)
xlabel('spectro parameters')
ylabel('Stims')
title('log(spectro)')
colorbar()
subplot(1,4,2);
imagesc(AvSpec_local)
xlabel('spectro parameters')
ylabel('Stims')
title('Average Spectro')
colorbar()
subplot(1,4,3);
imagesc(X_minusMean)
xlabel('spectro parameters')
ylabel('Stims')
title('Spectro - Mean(Spectro)')
colorbar()
subplot(1,4,4);
imagesc(x_ZC)
xlabel('spectro parameters')
ylabel('Stims')
title('Z-scored spectro')
colorbar()
pause(0.5)
            
    

