function [HalfTrain1, HalfTrain2, NumTrials]=organiz_data4coherence(Trials,Spectro,Win,ResDelay)
%% define new dataset depending on the size of the window of the model
% loop through the stims and only keep the Win first ms of them when
% they are longer than Win ms
NbStim=length(Trials);
HalfTrain1=[];
HalfTrain2=[];
NumTrials=nan(NbStim,1);
for dd = 1:NbStim
    NTrials=length(Trials{dd});
    duration=ceil(Spectro.to{dd}(end)*1000); %converting s in ms here
    if duration<=(Win+ResDelay)
        PSTH_local=zeros(NTrials,duration);
    else
        PSTH_local=zeros(NTrials,Win+ResDelay);
    end
    for tt=1:length(Trials{dd})
        Spike_times=Trials{dd}{tt}(find(Trials{dd}{tt}<= (Win + ResDelay)));
        ns = length(Spike_times);
        for nn=1:ns
            SpikeIdx=ceil(Spike_times(nn));
            if (SpikeIdx < 1 || SpikeIdx > duration)
                fprintf(1, 'Warning time index out of bounds for stim# %s trial %d: time_ind = %d ntimebins = %d\n', dd, tt, SpikeIdx, duration);
                continue;
            end
            PSTH_local(tt,SpikeIdx)=1+PSTH_local(tt,SpikeIdx);
        end
    end
    if mod(NTrials,2)==0
        Trialset1=1:2:NTrials;
        Trialset2=2:2:NTrials;
        NumTrials(dd)=NTrials;
    else
        Trialset1=1:2:(NTrials-1);
        Trialset2=2:2:(NTrials-1);
        NumTrials(dd)=NTrials-1;
    end
    HalfTrain1 = [HalfTrain1 mean(PSTH_local(Trialset1,:),1)];
    HalfTrain2 = [HalfTrain2 mean(PSTH_local(Trialset2,:),1)];
    
end

end