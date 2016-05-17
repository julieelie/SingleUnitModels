function [Info] = ModelInfoWrapper(ypredict)
Info = nan(length(ypredict),1);
for yy=1:length(ypredict)
    y_local = ypredict{yy};
    Info(yy) = info_model_Calculus(y_local);
end
end