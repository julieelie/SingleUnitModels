function [HalfTrain1, HalfTrain2]=organiz_SpikeArray4coherence(Spike_array,ParamModel)
%% define new dataset depending on the size of the Max window
% loop through the stims and only keep the Win first ms of the response when
% they are longer than MaxWin ms

% Identify how many sets of half sets of stims we want to obtain depending
% on the number of response duration that we want to test
Nb_MaxWin = length(ParamModel.MaxWin);

NbStim=length(Spike_array);
HalfTrain1=cell(Nb_MaxWin,1);
HalfTrain2=cell(Nb_MaxWin,1);
for MW = 1:Nb_MaxWin
    MaxWin_local = ParamModel.MaxWin(MW)* ParamModel.Response_samprate/1000;
    HalfTrain1{MW} = [];
    HalfTrain2{MW} = [];
    for dd = 1:NbStim
        NTrials = size(Spike_array{dd},1);
    
        if mod(NTrials,2)==0
            Trialset1=1:2:NTrials;
            Trialset2=2:2:NTrials;
        else
            Trialset1=1:2:(NTrials-1);
            Trialset2=2:2:(NTrials-1);
        end
        HalfTrain1{MW} = [HalfTrain1{MW} mean(Spike_array(Trialset1,1:MaxWin_local),1)];
        HalfTrain2{MW} = [HalfTrain2{MW} mean(Spike_array(Trialset2,1:MaxWin_local),1)];
    
    end
end
end