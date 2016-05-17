function [NeuroRes, PG_Index, FanoFactor_Index, Wins] = PoissonGaussianNeuralResponses(Spectro, VocType, PSTH, Trials,Cellname, MinWin, MaxWin, Increment, ResDelay, NeuroRes)
FIG=1;

if nargin<11
    NeuroRes = 'count';
end
if nargin<10
    ResDelay = 10; %predict the neural response with a 10ms delay after the end of the stimulus
end
if nargin<9
    Increment = 5; %increase the size of the spectro window with a 5ms pace
end
if nargin<8
    MaxWin = 150; %maximum values the window of analysis can reach
end
if nargin<7
    MinWin = 10; %minimum size of the window of analysis from the begining and also size of analysis of spike rate
end

%define the increasing size of the window of the spectrogram
Wins = MinWin:Increment:MaxWin;

% # of models to run on the data
modNum = length(Wins);

% Number of stims in the data set
NbStim = length(VocType);

%% Initialize output variables
PG_Index = nan(modNum,1);
FanoFactor_Index = nan(modNum,1);

%% Now loop through window sizes and look at spike rate distributions
for mm = 1:modNum
    %fprintf(1,'%d/%d models\n', mm, modNum);
    Win = Wins(mm);
    %% define new dataset depending on the size of the window of the model
    % loop through the stims and only keep the Win first ms of them when
    % they are longer than Win ms or disgard
    
    duration = nan(NbStim,1);
    for ss = 1:NbStim
        duration(ss)=Spectro.to{ss}(end)*1000; %converting s in ms here
    end
    Stim_local = find(duration >= (Win+ResDelay));% here we add ResDelay because we need to get sounds with corresponding psth that go ResDelay beyond the spectrogram of size Win
    NbStim_local = length(Stim_local);
    if NbStim_local<20
        sprintf('Only %d stims long enough to run the model: no model is run with window size %dms\n', NbStim_local, Win);
        break
    end
    %NBPC = [1:9 10:5:(NbStim_local*0.8)]; % This is not a good solution since the R2A profile of most cells show stairs-like structure with increasing number of PC we need to apply a ridge regression after the PCA.
    y = nan(NbStim_local,1);%this matrix will contain the average spike rate in spikes/ms at that precise position and for all the stims choosen for that run
    y_dev=cell(NbStim_local,1);
    ymean=nan(NbStim_local,1);
    yvar=ymean;
    for ss = 1:NbStim_local
        dd=Stim_local(ss);
        
        % Values of max spike rate and mean spike rate within the window
        % Check the distribution of responses (Gaussian or Poisson) for each stim
        if strcmp(NeuroRes, 'max')
            y(ss) = max(PSTH{dd}((Win-MinWin+ResDelay):(Win+ResDelay)));
        elseif strcmp(NeuroRes, 'mean')
            y(ss) = mean(PSTH{dd}((Win-MinWin+ResDelay):(Win+ResDelay)));% here we get the average Spike Rate over bins of 1ms so the spike rate is in spike/ms
        elseif strcmp(NeuroRes, 'count')
            y_dev{ss}=nan(length(Trials{dd}),1);
            for tt=1:length(Trials{dd})
                y_dev{ss}(tt)=sum((Trials{dd}{tt}>(Win-MinWin+ResDelay)).*(Trials{dd}{tt}<(Win+ResDelay)));
            end
            ymean(ss)=mean(y_dev{ss});
            yvar(ss)=var(y_dev{ss});
        else
            fprintf('please correctly write what kind of neural response you want to predict\n %s does not make any sense!!\n', NeuroRes);
    
        end
    end
    
    % Investigate how poisson or gaussian neural responses are for this
    % neuron
    MAX=max(max(yvar),max(ymean));
    PG_Index(mm) = sum(power(yvar-repmat(mean(yvar),length(yvar),1),2))/sum(power(yvar-ymean,2));
    FanoFactor_Index(mm) = nanmean(yvar./ymean);
    if FIG>0
        figure(1)
        plot(ymean,yvar,'.', 'MarkerSize',10)
        ylabel('Variance spike counts per stim')
        xlabel('mean spike count per stim')
        title(sprintf('Win=%d PoissonGaussian Index=%f\n FanoFactor=%f\n',Win,PG_Index(mm),FanoFactor_Index(mm)))
        hold on
        line([0 MAX], [0 MAX]);
        hold off
        pause()
    end    
         
end
if FIG>0
    figure()
    subplot(1,2,1)
    plot(1:modNum,log2(PG_Index))
    set(gca,'XTick',1:modNum)
    set(gca,'XTickLabel',Wins)
    ylabel('log2(PG_Index) >0 Poisson <0 Gaussian')
    xlabel('windows in ms')
    title(sprintf('SS Errors ratio Gaussian/Poison\n%s',Cellname));
    line([0 modNum], [0 0])
    subplot(1,2,2)
    plot(1:modNum,FanoFactor_Index)
    set(gca,'XTick',1:modNum)
    set(gca,'XTickLabel',Wins)
    ylabel('FanoFactor')
    xlabel('windows in ms')
    title(sprintf('Fano Factor var/mean\n%s',Cellname));
    pause(1)
end
end

