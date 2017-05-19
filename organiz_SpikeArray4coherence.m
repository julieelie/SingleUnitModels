function [HalfTrain1, HalfTrain2]=organiz_SpikeArray4coherence(Spike_array,ParamModel)
%% define new dataset depending on the size of the Max window
% loop through the stims and only keep the Win first ms of the response when
% they are longer than MaxWin ms


NbStim=length(Spike_array);

MaxWin_local = ParamModel.MaxWin* ParamModel.Response_samprate/1000;
HalfTrain1 = [];
HalfTrain2 = [];
for dd = 1:NbStim
    NTrials = size(Spike_array{dd},1);
    
    if mod(NTrials,2)==0
        Trialset1=1:2:NTrials;
        Trialset2=2:2:NTrials;
    else
        Trialset1=1:2:(NTrials-1);
        Trialset2=2:2:(NTrials-1);
    end
    HalfTrain1 = [HalfTrain1 mean(Spike_array{dd}(Trialset1,1:MaxWin_local),1)];
    HalfTrain2 = [HalfTrain2 mean(Spike_array{dd}(Trialset2,1:MaxWin_local),1)];
    
end
end
