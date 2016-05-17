function [Coefficients, CoefficientsNames, NumSignifdim, SSE, SST, R2] = myProgressiveLinearModel(y, x,Cat,RunAll,NSDim,varargin)
%This function progressively used the columns of x as parameters in the
%model that predict y with the categorical predictor Cat. Each column of x
%is used with Cat to fit sucessively the residuals of the previous model.
%If NSDim is not precised, all columns of x are considered and the
%algorithm stops using further column when the linear model is no more
%significant. If NSDim is precised the first NSDim columns of x are used in
%the model unless they are not significant.
%   y is a 1 X n vector of responses, x is a n X p1 matrix of predictors,
%   Cat is a 1 X n array of categorical predictor.
Interaction=0;
if nargin<4
    RunAll=0;
end

if nargin<5
    NSDim=size(x,2);
end

pnames = {'ValidatingSet'};
dflts  = {[]};
[ValSet] = internal.stats.parseArgs(pnames,dflts,varargin{:});
if isempty(ValSet)
    fprintf(1, 'No Validating Set: the sum of the squared errors will be calculated on the training dataset\n');
    ValSet.x=x;
    ValSet.y=y;
    ValSet.Cat=Cat;
end

Pval=0;
PvalAll = zeros(size(x,2),1);
ii=0;
Ucat = unique(Cat);
Ncat_local = length(Ucat);
Ind_Cat = 1:(Ncat_local+1);
Ind_Cat(2)=[];
Coeff.Cat=zeros(Ncat_local,1);
Coeff.X= zeros(size(x,2),1);
if Interaction==1
    Coeff.Int=zeros((Ncat_local-1)*size(x,2),1);
    CoeffNames.Int=cell((Ncat_local-1)*size(x,2),1);
end
CoeffNames.Cat=cell(Ncat_local,1);
CoeffNames.X=cell(size(x,2),1);

SSE = nan(size(x,2),1);
R2 = SSE;
y_local=y;
if RunAll
    DONE=0;
    while ii<size(x,2);
        ii=ii+1;
        ds3=dataset();
        ds3.(sprintf('DIM%d',ii)) = x(:,ii);
        ds3.Cattype=ordinal(Cat);
        ds3.y=y_local;
        if Interaction==1
            mdl3=fitlm(ds3,sprintf('y ~ DIM%d * Cattype',ii),'CategoricalVars',2);
        else
            mdl3=fitlm(ds3,sprintf('y ~ DIM%d + Cattype',ii),'CategoricalVars',2);
        end
        tbl3=anova(mdl3,'summary');
        Pval=tbl3.pValue(2);
        PvalAll(ii)=Pval;
        if Pval>0.05 && DONE==0
            DONE=1;
            NumSignifdim = ii-1;
        end
        y_local = y_local - mdl3.predict;
        Coeff.Cat = Coeff.Cat + mdl3.Coefficients.Estimate(Ind_Cat);
        CoeffNames.Cat = mdl3.Coefficients.Properties.ObsNames(Ind_Cat);
        Coeff.X(ii)=mdl3.Coefficients.Estimate(2);
        CoeffNames.X(ii) = mdl3.Coefficients.Properties.ObsNames(2);
        if Interaction==1
            Ind_Interaction = ((Ncat_local-1)*(ii-1)+1):((Ncat_local-1)*ii);
            Coeff.Int(Ind_Interaction)=mdl3.Coefficients.Estimate((Ind_Cat(end)+1):end);
            CoeffNames.Int(Ind_Interaction) = mdl3.Coefficients.Properties.ObsNames((Ind_Cat(end)+1):end);
        end
        
        % Calculate SSE for the training or validating data set
        y_predicted_X = ValSet.x(:,1:ii) * Coeff.X(1:ii,1);
        y_predicted_Cat = nan(size(ValSet.x,1),1);
        if Interaction==1
            y_predicted_Int = nan(size(ValSet.x,1),1);
            for ts=1:size(ValSet.x,1)
                Ind_CatType = strcmp(CoeffNames.Cat,sprintf('Cattype_%s',ValSet.Cat{ts}));
                Ind_int = zeros(size(CoeffNames.Int));
                for nn = 1:ii
                    Ind_int = Ind_int + strcmp(CoeffNames.Int, sprintf('DIM%d:Cattype_%s',nn,ValSet.Cat{ts}));
                end
                y_predicted_Cat(ts) = Coeff.Cat(1) + sum(Coeff.Cat.*Ind_CatType);
                if sum(Ind_int)>0
                    y_predicted_Int(ts) = ValSet.x(ts,1:ii)*Coeff.Int(find(Ind_int));
                else
                    y_predicted_Int(ts) = 0;
                end
            end
            y_predicted = y_predicted_X + y_predicted_Cat + y_predicted_Int;
        else
            for ts=1:size(ValSet.x,1)
                Ind_CatType = strcmp(CoeffNames.Cat,sprintf('Cattype_%s',ValSet.Cat{ts}));
                y_predicted_Cat(ts) = Coeff.Cat(1) + sum(Coeff.Cat.*Ind_CatType);
            end
            y_predicted = y_predicted_X + y_predicted_Cat;
        end
            
        
        SSE(ii) = sum((ValSet.y - y_predicted).^2);
        SST_local = sum((ValSet.y-mean(y)).^2);
        R2(ii) = 1 - SSE(ii)./SST_local;
    end
    SST=repmat(SST_local, size(x,2),1);
    figure(9)
    subplot(1,2,1)
    plot(1:ii,SSE)
    ylabel('SSE');
    xlabel('Number of Dimensions of x used');
    subplot(1,2,2)
    plot(1:ii,R2);
    ylabel('Rsquare on same data set');
    xlabel('Number of Dimensions of x used');
else %The progressive model stop incorporating new dimensions of x when the model is no longer significant
    while Pval<=0.05 && ii<NSDim
        ii=ii+1;
        ds3=dataset();
        ds3.(sprintf('DIM%d',ii)) = x(:,ii);
        ds3.Cattype=ordinal(Cat);
        ds3.y=y_local;
        if Interaction==1
            mdl3=fitlm(ds3,sprintf('y ~ DIM%d * Cattype',ii),'CategoricalVars',2);
        else
            mdl3=fitlm(ds3,sprintf('y ~ DIM%d + Cattype',ii),'CategoricalVars',2);
        end
        tbl3=anova(mdl3,'summary');
        Pval=tbl3.pValue(2);
        if Pval>0.05
            NumSignifdim = ii-1;
            break
        else
            NumSignifdim = NSDim;
        end
        y_local = y_local - mdl3.predict;
        Coeff.Cat = Coeff.Cat + mdl3.Coefficients.Estimate(Ind_Cat);
        CoeffNames.Cat = mdl3.Coefficients.Properties.ObsNames(Ind_Cat);
        Coeff.X(ii)=mdl3.Coefficients.Estimate(2);
        CoeffNames.X(ii) = mdl3.Coefficients.Properties.ObsNames(2);
        if Interaction==1
            Ind_Interaction = ((Ncat_local-1)*(ii-1)+1):((Ncat_local-1)*ii);
            Coeff.Int(Ind_Interaction)=mdl3.Coefficients.Estimate((Ind_Cat(end)+1):end);
            CoeffNames.Int(Ind_Interaction) = mdl3.Coefficients.Properties.ObsNames((Ind_Cat(end)+1):end);
        end
        % Calculate SSE for the training or validating data set
        y_predicted_X = ValSet.x(:,1:ii) * Coeff.X(1:ii,1);
        y_predicted_Cat = nan(size(ValSet.x,1),1);
        if Interaction==1
            y_predicted_Int = nan(size(ValSet.x,1),1);
            for ts=1:size(ValSet.x,1)
                Ind_CatType = strcmp(CoeffNames.Cat,sprintf('Cattype_%s',ValSet.Cat{ts}));
                Ind_int = zeros(size(CoeffNames.Int));
                for nn = 1:ii
                    Ind_int = Ind_int + strcmp(CoeffNames.Int, sprintf('DIM%d:Cattype_%s',nn,ValSet.Cat{ts}));
                end
                y_predicted_Cat(ts) = Coeff.Cat(1) + sum(Coeff.Cat.*Ind_CatType);
                if sum(Ind_int)>0
                    y_predicted_Int(ts) = ValSet.x(ts,1:ii)*Coeff.Int(find(Ind_int));
                else
                    y_predicted_Int(ts) = 0;
                end
            end
            y_predicted = y_predicted_X + y_predicted_Cat + y_predicted_Int;
        else
            for ts=1:size(ValSet.x,1)
                Ind_CatType = strcmp(CoeffNames.Cat,sprintf('Cattype_%s',ValSet.Cat{ts}));
                y_predicted_Cat(ts) = Coeff.Cat(1) + sum(Coeff.Cat.*Ind_CatType);
            end
            y_predicted = y_predicted_X + y_predicted_Cat;
        end
        SSE(ii) = sum((ValSet.y - y_predicted).^2);
        SST_local = sum((ValSet.y-mean(y)).^2);
        R2(ii) = 1 - SSE(ii)./SST_local;
    end
    SSE=SSE(1:NumSignifdim);
    R2=R2(1:NumSignifdim);
end
if Interaction==1
    Coefficients=[Coeff.Cat;Coeff.X(1:NumSignifdim);Coeff.Int(1:((Ncat_local-1)*NumSignifdim))]';
    CoefficientsNames={CoeffNames.Cat{:} CoeffNames.X{1:NumSignifdim} CoeffNames.Int{1:((Ncat_local-1)*NumSignifdim)}};
else
    Coefficients=[Coeff.Cat;Coeff.X(1:NumSignifdim)]';
    CoefficientsNames={CoeffNames.Cat{:} CoeffNames.X{1:NumSignifdim}};
end
end

