function [ValSet, TestSet] = create_cross_validation_sets(VocType, Emittername)
%% Construct a testing dataset and a validating dataset
% remove one emitter per call category, if one category contain only
% one emitter, don't use it in the model
CT = unique(VocType);
ValSet=nan(size(VocType));
TestSet=nan(size(VocType));
jj=1;
kk=1;
for Cc=1:length(CT);
    Ccs = CT(Cc);
    Indcc = find(strcmp(VocType,Ccs));
    Em = unique(Emittername(Indcc));
    if length(Em)>1
        RandInd = randperm(length(Em));
        if ~strcmp(Em(RandInd(1)),'STRFxx0000')
            IndEm = find(strcmp(Emittername, Em(RandInd(1))));
        else
            IndEm = find(strcmp(Emittername, Em(RandInd(2))));
        end
        IndccEmVal = intersect(IndEm, Indcc);
        IndccEmTes = setdiff(Indcc, IndccEmVal);
        ValSet(jj:(jj+length(IndccEmVal)-1))=IndccEmVal;
        jj = jj+length(IndccEmVal);
        TestSet(kk:(kk+length(IndccEmTes)-1))=IndccEmTes;
        kk = kk+length(IndccEmTes);
    else
        fprintf(1, '%s is not used to calculate the model (only one emitter)', Ccs);
    end
end
ValSet=ValSet(1:(jj-1));
TestSet = TestSet(1:(kk-1));

end

