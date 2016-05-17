STRF_show=0;
NL_show=0;
LL_show=0;
%load('SynCell_CellResultsAlphaLambdaDev_AllModels_poisson_log.mat')
%load('DC128_GLMPoisson.mat')
%load('Ag127_CellResultsAlphaLambdaDev_AllModels_poisson_log.mat')
%load('ACCell_CellResultsAlphaLambdaDev_AllModels_poisson_log.mat')
% load('WholeVoc_Site2_L1000R900_e13_s0_ss1_GLMPoisson.mat')% SYN
% load('WholeVoc_Site2_L1100R1450_e14_s0_ss1_GLMPoisson.mat')% Ac
% load('WholeVoc_Site3_L2500R2300_e22_s1_ss1_GLMPoisson.mat') %Ag127
% load('WholeVoc_Site4_L1500R1900_e23_s0_ss2_GLMPoisson.mat')%DC128
cd /Users/elie/Documents/CODE/data/matfile/ModMat
PoissonFiles=dir('Models*');

MinDevlog10Lambdas.Acoustic = [];
MinDevlog10Lambdas.AcSem = [];
MinDevlog10Lambdas.AcSem2 = [];

for ff=1:length(PoissonFiles)
    load(PoissonFiles(ff).name)
    %% Loop through structure and plot boxplot of best values of Deviance for each alphas
    Models = fieldnames(Deviance);
    for mm=1:modNum
        for ii=1:numel(Models)
            if strcmp(Models{ii}(1),'A') % Make sure this is a Acoustic or Acsem model
                Deviance_All = Deviance.(Models{ii});

                % Find Max min values of Deviance in the dataset
                BT=length(Deviance_All.lambda{mm});
                NAlphas = size(Deviance_All.lambda,2);
                MaxDev=nan(NAlphas*BT,1);
                MinDev=nan(NAlphas*BT,1);
                LambdaMinDev = nan(NAlphas*BT,1);
                Lambda_values=nan(NAlphas*BT*100,1);
                Deviance_values = Lambda_values;
                Alpha_values = Lambda_values;
                dd=0;
                vv=0;
                for AA=1:NAlphas
                    Deviance_local=Deviance_All.values{mm,AA};
                    Lambda_local=Deviance_All.lambda{mm,AA};
                    for BB=1:length(Deviance_local)
                        dd=dd+1;
                        MaxDev(dd)=max(Deviance_local{BB});
                        MinDev(dd)=min(Deviance_local{BB});
                        LambdaMinDev(dd)=log10(Lambda_local{BB}(Deviance_local{BB}==min(Deviance_local{BB})));
                        for ll=1:length(Deviance_local{BB})
                            vv=vv+1;
                            Deviance_values(vv)=Deviance_local{BB}(ll);
                            Lambda_values(vv)=log10(Lambda_local{BB}(ll));
                            Alpha_values(vv)=Alphas(AA);
                        end
                    end
                end
                MinDev=MinDev(1:dd);
                LambdaMinDev=LambdaMinDev(1:dd);
                MAXDev = max(MaxDev);
                MINDev = min(MinDev);
                Alpha_values=Alpha_values(1:vv);
                Deviance_values=Deviance_values(1:vv);
                Lambda_values=Lambda_values(1:vv);
                %cubehelix_niceplot(Alpha_values,Lambda_values,Deviance_values)
        %         figure(ii)
        %         plot3(Alpha_values,Lambda_values,Deviance_values,'.')
        %         xlabel('Alpha')
        %         ylabel('Log lambda')
        %         zlabel('Deviance')
        %         set(gca,'ZLim',[mean(Deviance_values)-2.*std(Deviance_values) mean(Deviance_values)+3.*std(Deviance_values)/(length(Deviance_values))^0.5])
        %         %set(gca,'ZLim',[0 1])
        %         title(sprintf('%d bootstraps %d lambda values %s model',BT,length(Deviance_local),Models{ii}))
                % Noise in lambda minimum values and deviance minimumu values
        %         figure(ii)
        %         boxplot(Deviance_All.bestvalues{mm})
        %         xlabel('Alpha')
        %         ylabel('Minimum Deviance')
        %         set(gca,'XTickLabel',Alphas)
        %         title(sprintf('%s model',Models{ii}));
            end
        end

    % figure(5)
    % boxplot(Deviance.Sem.bestvalues{mm})
    % xlabel('Alpha')
    % ylabel('Minimum Deviance')
    % set(gca,'XTickLabel',Alphas)
    % title('Semantic model');

        figure(4)
        boxplot([Deviance.Acoustic.bestvalues{mm} Deviance.AcSem2.bestvalues{mm} Deviance.AcSem.bestvalues{mm} Deviance.Sem.bestvalues{mm}])
        xlabel(sprintf('Models win=%d',Wins(mm)))
        ylabel('Minimum Deviance')
        set(gca,'XTickLabel',{sprintf('Acoustic FitIndex=%f',Deviance.Acoustic.FitIndex{mm}),sprintf('AcSem AcOffset FitIndex=%f',Deviance.AcSem2.FitIndex{mm}),sprintf('AcSem SemOffset FitIndex=%f',Deviance.AcSem.FitIndex{mm}),sprintf('Semantic FitIndex=%f',Deviance.Sem.FitIndex{mm})})
        pause()
    end

    %% Ploting loglikelihood as a function of windows if we have more than one window
    if LL_show
        if modNum>1 && NAlphas==1
            LLAverage.Acoustic = nan(modNum,1);
            LLAverage.AcSem = nan(modNum,1);
            LLAverage.AcSem2 = nan(modNum,1);
            LLAverage.Sem = nan(modNum,1);
            LLAverage.Ceiling = nan(modNum,1);
            LLAverage.Floor = nan(modNum,1);
            LLAverage.AutoRegressive = nan(modNum,1);
            for mm=1:modNum
                LLAverage.Acoustic(mm) = mean(LL.Acoustic.bestvalues{mm,1});
                LLAverage.AcSem(mm) =mean(LL.AcSem.bestvalues{mm,1});
                LLAverage.AcSem2(mm) = mean(LL.AcSem2.bestvalues{mm,1});
                LLAverage.Sem(mm) = mean(LL.Sem.bestvalues{mm});
                LLAverage.Ceiling(mm) = mean(LL.Ceiling.bestvalues{mm,1});
                LLAverage.Floor(mm) = mean(LL.Floor.bestvalues{mm,1});
                LLAverage.AutoRegressive(mm) = mean(LL.AutoRegressive.bestvalues{mm,1});
            end
            models=fieldnames(LLAverage);
            Colors='brmgcky';
            figure(5)
            for ii=1:numel(models)
                LLAverage_local=LLAverage.(models{ii});
                plot(LLAverage_local,sprintf(Colors(ii)))
                hold on
            end
            legend(models);
            ylabel('LogLikelihood');
            xlabel('Window position');
            set(gca,'XTick',[1 2 3 4],'XTickLabel', Wins,'XLim',[0 modNum+1]);
            hold off

        end
    end
    %% Testing differences between models based on Deviance calculated on validating dataset
    for mm=1:modNum
        % Acoustic vs Acsem1 Model
        DiffDevianceAcousticAcSem=(Deviance.Acoustic.bestvalues{mm} - Deviance.AcSem.bestvalues{mm})./PropVal.values{mm};
        %one sample t-test
        [h,pAcousticAcSem,ci,stats]=ttest(DiffDevianceAcousticAcSem);

        % Semantic vs Acsem1 Model
        DiffDevianceSemanticAcSem=(Deviance.Sem.bestvalues{mm} - Deviance.AcSem.bestvalues{mm})./PropVal.values{mm};
        %one sample t-test
        [h,pSemanticAcSem,ci,stats]=ttest(DiffDevianceSemanticAcSem);

        % Semantic vs Acoustic Model
        DiffDevianceAcousticSemantic=(Deviance.Acoustic.bestvalues{mm} - Deviance.Sem.bestvalues{mm})./PropVal.values{mm};
        %one sample t-test
        [h,pAcousticSemantic,ci,stats]=ttest(DiffDevianceAcousticSemantic);

        % Acoustic vs Acsem2 Model
        DiffDevianceAcousticAcSem2=(Deviance.Acoustic.bestvalues{mm} - Deviance.AcSem2.bestvalues{mm})./PropVal.values{mm};
        %one sample t-test
        [h,pAcousticAcSem2,ci,stats]=ttest(DiffDevianceAcousticAcSem2);

        % Semantic vs Acsem2 Model
        DiffDevianceSemanticAcSem2=(Deviance.Sem.bestvalues{mm} - Deviance.AcSem2.bestvalues{mm})./PropVal.values{mm};
        %one sample t-test
        [h,pSemanticAcSem2,ci,stats]=ttest(DiffDevianceSemanticAcSem2);

    figure(6)
    subplot(2,3,1)
    boxplot(DiffDevianceAcousticAcSem)
    xlabel('Alpha')
    ylabel('Difference of Deviance (Best guess Model)')
    set(gca,'XTickLabel',Alphas)
    title(sprintf('Acoustic - AcSem SemOffset\nttest pvalue=%f',pAcousticAcSem));
    subplot(2,3,2)
    boxplot(DiffDevianceSemanticAcSem)
    xlabel('Alpha')
    ylabel('Difference of Deviance (Best guess Model)')
    set(gca,'XTickLabel',Alphas)
    title(sprintf('Semantic - AcSem SemOffset\nttest pvalue=%f',pSemanticAcSem));
    subplot(2,3,3)
    boxplot(DiffDevianceAcousticSemantic)
    xlabel('Alpha')
    ylabel('Difference of Deviance (Best guess Model)')
    set(gca,'XTickLabel',Alphas)
    title(sprintf('Acoustic - Semantic\nttest pvalue=%f',pAcousticSemantic));
    subplot(2,3,4)
    boxplot(DiffDevianceAcousticAcSem2)
    xlabel('Alpha')
    ylabel('Difference of Deviance (Best guess Model)')
    set(gca,'XTickLabel',Alphas)
    title(sprintf('Acoustic - AcSem Ac Offset\nttest pvalue=%f',pAcousticAcSem2));
    subplot(2,3,5)
    boxplot(DiffDevianceSemanticAcSem2)
    xlabel('Alpha')
    ylabel('Difference of Deviance (Best guess Model)')
    set(gca,'XTickLabel',Alphas)
    title(sprintf('Semantic - AcSem AcOffset\nttest pvalue=%f',pSemanticAcSem2));
    pause()
    end
    % subplot(2,3,4)
    % boxplot(-DiffAIC_AcSemAcoustic)
    % xlabel('Alpha')
    % ylabel('Difference of AIC')
    % set(gca,'XTickLabel',Alphas)
    % title('Acoustic - AcSem');
    % subplot(2,3,5)
    % boxplot(-DiffAIC_AcSemSemantic)
    % xlabel('Alpha')
    % ylabel('Difference of AIC')
    % set(gca,'XTickLabel',Alphas)
    % title('Semantic- AcSem');
    
    
    %% Find the best value of lambdas over bootstrap by summing deviance
    SumDev.Acoustic = zeros(modNum,BT);
    SumDev.AcSem = SumDev.Acoustic;
    SumDev.AcSem2 = SumDev.Acoustic;
    for mm=1:modNum
        AA=1;
        MAXL_Ac=[];
        MINL_Ac=[];
        MAXL_AcSem=[];
        MINL_AcSem=[];
        MAXL_AcSem2=[];
        MINL_AcSem2=[];
        for BB=1:BT
            MAXL_Ac=max([MAXL_Ac Deviance.Acoustic.lambda{mm,AA}{BB}]);
            MINL_Ac=min([MINL_Ac Deviance.Acoustic.lambda{mm,AA}{BB}]);
            MAXL_AcSem=max([MAXL_AcSem Deviance.AcSem.lambda{mm,AA}{BB}]);
            MINL_AcSem=min([MINL_AcSem Deviance.AcSem.lambda{mm,AA}{BB}]);
            MAXL_AcSem2=max([MAXL_AcSem2 Deviance.AcSem2.lambda{mm,AA}{BB}]);
            MINL_AcSem2=min([MINL_AcSem2 Deviance.AcSem2.lambda{mm,AA}{BB}]);
        end
        L_Ac = linspace(log10(MINL_Ac),log10(MAXL_Ac),BT);
        L_AcSem = linspace(log10(MINL_AcSem),log10(MAXL_AcSem),BT);
        L_AcSem2 = linspace(log10(MINL_AcSem2),log10(MAXL_AcSem2),BT);



        for BB=1:BT
            vq = interp1(log10(Deviance.Acoustic.lambda{mm,AA}{BB}),Deviance.Acoustic.values{mm,AA}{BB}/PropVal.values{mm}(BB,AA),L_Ac);
            SumDev.Acoustic(mm,:) = SumDev.Acoustic(mm,:) + vq;
            vq = interp1(log10(Deviance.AcSem.lambda{mm,AA}{BB}),Deviance.AcSem.values{mm,AA}{BB}/PropVal.values{mm}(BB,AA),L_AcSem);
            SumDev.AcSem(mm,:) = SumDev.AcSem(mm,:) + vq;
            vq = interp1(log10(Deviance.AcSem2.lambda{mm,AA}{BB}),Deviance.AcSem2.values{mm,AA}{BB}/PropVal.values{mm}(BB,AA),L_AcSem2);
            SumDev.AcSem2(mm,:) = SumDev.AcSem2(mm,:) + vq;
        end

        MinDevlog10Lambdas.Acoustic = [MinDevlog10Lambdas.Acoustic L_Ac(find(SumDev.Acoustic(mm,:)==min(SumDev.Acoustic(mm,:))))];
        MinDevlog10Lambdas.AcSem = [MinDevlog10Lambdas.AcSem L_AcSem(find(SumDev.AcSem(mm,:)==min(SumDev.AcSem(mm,:))))];
        MinDevlog10Lambdas.AcSem2 = [MinDevlog10Lambdas.AcSem2 L_AcSem2(find(SumDev.AcSem2(mm,:)==min(SumDev.AcSem2(mm,:))))];

        figure(7)
        plot(L_Ac,SumDev.Acoustic(mm,:),'k*-',L_AcSem,SumDev.AcSem(mm,:),'r*-',L_AcSem,SumDev.AcSem2(mm,:),'y*-')
        legend('Acoustic model','AcSem Sem Offset Model','AcSem Ac Offset Model')
        ylabel(sprintf('sum of deviance per obs over %d bootstrap',BT))
        xlabel('log10 lambda')
        title(sprintf('%s Win=%d',Cellname, Wins(mm)))
        hold on
        pause()
    end
    hold off


    %% Plot best STRFs for each alpha each model, each boostrap
    % Also calculate an R2
    R2.SS_Acoustic=nan(BootstrapSTRF,NAlphas);
    R2.SS_AcSem = R2.SS_Acoustic;
    R2.SS_Semantic = R2.SS_Acoustic;
    R2.SS_T = R2.SS_Acoustic;
    Av_Lambda.Acoustic = nan(BootstrapSTRF,NAlphas);
    Av_Lambda.AcSem = Av_Lambda.Acoustic;
    Av_Deviance.Acoustic = nan(BootstrapSTRF,NAlphas);
    Av_Deviance.AcSem=Av_Deviance.Acoustic;
    MinLambdas.AcSem=nan(BootstrapSTRF,NAlphas);
    MinLambdas.Acoustic=nan(BootstrapSTRF,NAlphas);
    
    
    for AA=1:1
        for BB=1:BT
            if STRF_show==1
                % plot Acoustic STRF
                figure(10)
                valc_Ac = max(abs(max(max(ModelB_Ac{BB,AA}))), abs(min(min(ModelB_Ac{BB,AA}))));
                if valc_Ac==0
                    imagesc(Model.TickSTRFspectro.to{mm}, Model.TickSTRFspectro.to{mm},ModelB_Ac{BB,AA})
                else
                    imagesc(Model.TickSTRFspectro.to{mm}, Model.TickSTRFspectro.to{mm},ModelB_Ac{BB,AA}, [-valc_Ac valc_Ac])
                end
                axis xy
                title(sprintf('Acoustic Alpha=%f %d/%d\nDeviance=%f Biais=%f',Alphas(AA),BB,BT,Deviance.Acoustic.bestvalues{mm}(BB,AA),ModelB0_Ac{BB,AA}));

                % Plot Model coefficients AcSem
                figure(11)
                subplot(1,2,1)
                valc_AcSem = max(abs(max(max(ModelBspectro_AcSem{BB,AA}))), abs(min(min(ModelBspectro_AcSem{BB,AA}))));
                if valc_AcSem==0
                    imagesc(Model.TickSTRFspectro.to{mm}, Model.TickSTRFspectro.to{mm},ModelBspectro_AcSem{BB,AA})
                else
                    imagesc(Model.TickSTRFspectro.to{mm}, Model.TickSTRFspectro.to{mm},ModelBspectro_AcSem{BB,AA}, [-valc_AcSem valc_AcSem])
                end
                axis xy
                title(sprintf('AcSem Alpha=%f %d/%d\nDeviance=%f Biais=%f',Alphas(AA),BB,BT,Deviance.AcSem.bestvalues{mm}(BB,AA),ModelB0_AcSem{BB,AA}));
                ss=subplot(1,2,2);
                plot([ModelB0_AcSem{BB,AA}; ModelBsem_AcSem{BB,AA}] + [0; repmat(ModelB0_AcSem{BB,AA},length(ModelBsem_AcSem{BB,AA}),1)],1:length(UVOC));
                title('Coefficients Categories');
                set(ss,'YTickLabel', UVOC);

                % Plot model coefficients Semantic
                figure(12)
                plot(ModelB_Sem{BB,AA} + [0 ; repmat(ModelB_Sem{BB,AA}(1),length(ModelB_Sem{BB,AA})-1,1)],1:length(UVOC));
                title(sprintf('Sem Model\ndeviance=%f', Deviance.Sem.bestvalues{mm}(BB,AA)));
                set(gca,'YTickLabel', UVOC);

                % Plot prediction for each model
                LocalModelPredict=cell(3,1);
                LocalModelPredict{2}=Model.AcSem.ypredictVal{mm,AA}{BB};
                LocalModelPredict{1}=Model.Acoustic.ypredictVal{mm,AA}{BB};
                LocalModelPredict{3}=Model.Semantic.ypredictVal{mm,AA}{BB};
                R2.SS_T(BB,AA) = sum(power(Model.yobsVal{mm,AA}{BB}-repmat(mean(Model.yobsVal{mm,AA}{BB}),length(Model.yobsVal{mm,AA}{BB}),1),2));
                R2.SS_Acoustic(BB,AA) = sum(power(Model.yobsVal{mm,AA}{BB}-Model.Acoustic.ypredictVal{mm,AA}{BB},2));
                R2.SS_Semantic(BB,AA) = sum(power(Model.yobsVal{mm,AA}{BB}-Model.Semantic.ypredictVal{mm,AA}{BB},2));
                R2.SS_AcSem(BB,AA) = sum(power(Model.yobsVal{mm,AA}{BB}-Model.AcSem.ypredictVal{mm,AA}{BB},2));
                R2.Acoustic(BB,AA)=1-R2.SS_Acoustic(BB,AA)/R2.SS_T(BB,AA);
                R2.Semantic(BB,AA)=1-R2.SS_Semantic(BB,AA)/R2.SS_T(BB,AA);
                R2.AcSem(BB,AA)=1-R2.SS_AcSem(BB,AA)/R2.SS_T(BB,AA);
                Legend= cell(3,1);
                Legend{1}=sprintf('Acoustic only R2=%f',1-R2.SS_Acoustic(BB,AA)/R2.SS_T(BB,AA));
                Legend{3}=sprintf('Semantic only R2=%f',1-R2.SS_Semantic(BB,AA)/R2.SS_T(BB,AA));
                Legend{2}=sprintf('Acoustic + Semantic R2=%f', 1-R2.SS_AcSem(BB,AA)/R2.SS_T(BB,AA));
                MAXX=max(Model.yobsVal{AA}{BB});
                MAXY=max([max(LocalModelPredict{1}) max(LocalModelPredict{2}) max(LocalModelPredict{3})]);
                MAX=max(MAXX,MAXY);

                for jj=1:3
                    figure(13)
                    subplot(1,3,jj);
                    h=gscatter(Model.yobsVal{mm,AA}{BB}, LocalModelPredict{jj}, Model.yobsCat{mm,AA}{BB}, 'mgcbrkyyr', '......d.d',[20 20 20 20 20 20 10 20 10]);
                    title(sprintf('%s', Legend{jj}));
                    xlabel('observed spike rate (spike per ms)')
                    ylabel('predicted spike rate (spike per ms)')
                    axis([0 MAX+MAX/4 0 MAX]);
                    hold on
                    plot(0:MAX/10:MAX,0:MAX/10:MAX, 'k');
                    hold off
                end
            end

            % calculate the average best lambda and Deviance over boostrap
            Av_Deviance.Acoustic(BB,AA)=mean(Deviance.Acoustic.bestvalues{mm}(1:BB,AA));
            MinLambdas.Acoustic(BB,AA)=Deviance.Acoustic.lambda{mm,AA}{BB}(find(Deviance.Acoustic.values{mm,AA}{BB}==Deviance.Acoustic.bestvalues{mm}(BB,AA)));
            Av_Lambda.Acoustic(BB,AA)=median(MinLambdas.Acoustic(1:BB,AA));
            Av_Deviance.AcSem(BB,AA)=mean(Deviance.AcSem.bestvalues{mm}(1:BB,AA));
            MinLambdas.AcSem(BB,AA)=Deviance.AcSem.lambda{mm,AA}{BB}(find(Deviance.AcSem.values{mm,AA}{BB}==Deviance.AcSem.bestvalues{mm}(BB,AA)));
            Av_Lambda.AcSem(BB,AA)=median(MinLambdas.AcSem(1:BB,AA));
            
            if STRF_show==1
                figure(9)
                plot(log10(Deviance.Acoustic.lambda{mm,AA}{BB}),Deviance.Acoustic.values{mm,AA}{BB},'k-',log10(Deviance.AcSem.lambda{mm,AA}{BB}),Deviance.AcSem.values{mm,AA}{BB},'r-');
                hold on
                plot(log10(MinLambdas.Acoustic(BB,AA)),Deviance.Acoustic.bestvalues{mm}(BB,AA),'k*')
                hold on
                plot(log10(MinLambdas.AcSem(BB,AA)),Deviance.AcSem.bestvalues{mm}(BB,AA),'r*')
                hold on
                hline(Deviance.Sem.bestvalues{mm}(BB,AA),'g-')
                legend('Acoustic', 'AcSem', 'Min Dev Acoustic','Min Dev AcSem')
                xlabel('log10 lambda')
                ylabel('Deviance')
                hold off
                pause(1)
            end
        end
    end
    %% Plot the cumulative average best Deviance and best lambda along bootstrap
    Color='mgcbrk';

    figure(14)
    for AA=1:length(Alphas)
        subplot(2,2,1)
        hold on
        plot(1:BootstrapSTRF,Av_Deviance.Acoustic(:,AA),sprintf('%s-',Color(AA)));
        subplot(2,2,2)
        hold on
        plot(1:BootstrapSTRF,Av_Deviance.AcSem(:,AA),sprintf('%s-',Color(AA)));
        subplot(2,2,3)
        hold on
        plot(1:BootstrapSTRF,log10(Av_Lambda.Acoustic(:,AA)),sprintf('%s-',Color(AA)));
        subplot(2,2,4)
        hold on
        plot(1:BootstrapSTRF,log10(Av_Lambda.AcSem(:,AA)),sprintf('%s-',Color(AA)));
    end
    subplot(2,2,1)
    ylabel('Average Deviance Acoustic Model')
    xlabel('Booststrap')
    subplot(2,2,2)
    ylabel('Average Deviance AcSem Model')
    xlabel('Booststrap')
    subplot(2,2,3)
    ylabel('Median Lambda Acoustic Model (log10)')
    xlabel('Booststrap')
    subplot(2,2,4)
    ylabel('Median Lambda AcSem Model (log10)')
    xlabel('Booststrap')
    hold off
    figure(14)
    pause()
    
    if NAlphas>1
        %% Plot deviance of Ac vs Acsem + diag and color per alpha see where it's
        % most often above/below diag chose alpha
        Color='mgcbrk';
        extremeBestDeviance=nan(NAlphas,2);
        legendcolor=cell(NAlphas,2);
        H=nan(NAlphas,2);
        for AA=1:NAlphas
            figure(15)
            subplot(1,2,1)
            hold on
            H(AA,1)=plot(Deviance.AcSem.bestvalues{mm}(:,AA),Deviance.Acoustic.bestvalues{mm}(:,AA),sprintf('%s.',Color(AA)), 'MarkerSize', 20);
            extremeBestDeviance(AA,1)=min(min(Deviance.Acoustic.bestvalues{mm}(:,AA)),min(Deviance.AcSem.bestvalues{mm}(:,AA)));
            extremeBestDeviance(AA,2)=max(max(Deviance.Acoustic.bestvalues{mm}(:,AA)),max(Deviance.AcSem.bestvalues{mm}(:,AA)));
            legendcolor{AA,1}=num2str(Alphas(AA));
            subplot(1,2,2)
            hold on
            %H(AA,2)=plot(mean(Deviance.AcSem.bestvalues{1}(:,AA)),mean(Deviance.Acoustic.bestvalues{1}(:,AA)),sprintf('%s.',Color(AA)), 'MarkerSize', 20);
            X=[mean(Deviance.AcSem.bestvalues{mm}(:,AA))-std(Deviance.AcSem.bestvalues{mm}(:,AA)) mean(Deviance.AcSem.bestvalues{mm}(:,AA))+std(Deviance.AcSem.bestvalues{mm}(:,AA))];
            Y=[mean(Deviance.Acoustic.bestvalues{mm}(:,AA))-std(Deviance.Acoustic.bestvalues{mm}(:,AA)) mean(Deviance.Acoustic.bestvalues{mm}(:,AA))+std(Deviance.Acoustic.bestvalues{mm}(:,AA))];
            H(AA,2)=line([mean(Deviance.AcSem.bestvalues{mm}(:,AA)) mean(Deviance.AcSem.bestvalues{mm}(:,AA))],Y,'Color',Color(AA));
            H(AA,2)=line(X,[mean(Deviance.Acoustic.bestvalues{mm}(:,AA)) mean(Deviance.Acoustic.bestvalues{mm}(:,AA))],'Color',Color(AA));
            legendcolor{AA,2}=sprintf('mean %s',num2str(Alphas(AA)));
        end
        subplot(1,2,1)
        legend(H(:,1),legendcolor{:,1})
        MAX = max(extremeBestDeviance(:,2));
        MIN = min(extremeBestDeviance(:,1));
        axis([MIN MAX+0.1 MIN MAX+0.1])
        hold on
        plot(MIN:MAX/10:(MAX+0.1),MIN:MAX/10:(MAX+0.1), 'k');
        hold off
        ylabel('Acoustic Model Deviance')
        xlabel('AcSem Model Deviance')
        subplot(1,2,2)
        legend(H(:,2),legendcolor{:,2})
        MAX = max(extremeBestDeviance(:,2));
        MIN = min(extremeBestDeviance(:,1));
        axis([MIN MAX+0.1 MIN MAX+0.1])
        hold on
        plot(MIN:MAX/10:(MAX+0.1),MIN:MAX/10:(MAX+0.1), 'k');
        hold off
        ylabel('Acoustic Model Deviance')
        xlabel('AcSem Model Deviance')
    end
    
    %% Correlation between deviance (AIC) and an r2??
    if NAlphas>1
        Color='mgcbrk';
        extremeBestDeviance=nan(NAlphas,2);
        legendcolor=cell(NAlphas,2);
        H=nan(NAlphas,3);
        for AA=1:NAlphas
            figure(16)
            subplot(2,2,1)
            hold on
            H(AA,1)=plot(Deviance.AcSem.bestvalues{mm}(:,AA),1-R2.SS_AcSem(:,AA)./R2.SS_T(:,AA),sprintf('%s.',Color(AA)), 'MarkerSize', 20);
            %extremeBestDeviance(AA,1)=min(min(Deviance.AcSem.bestvalues{1}(:,AA)),min(Deviance.AcSem.bestvalues{1}(:,AA)));
            %extremeBestDeviance(AA,2)=max(max(Deviance.AcSem.bestvalues{1}(:,AA)),max(Deviance.AcSem.bestvalues{1}(:,AA)));
            legendcolor{AA,1}=num2str(Alphas(AA));
            subplot(2,2,2)
            hold on
            H(AA,2)=plot(Deviance.Acoustic.bestvalues{mm}(:,AA),1-R2.SS_Acoustic(:,AA)./R2.SS_T(:,AA),sprintf('%s.',Color(AA)), 'MarkerSize', 20);
            %extremeBestDeviance(AA,1)=min(min(Deviance.Acoustic.bestvalues{1}(:,AA)),min(Deviance.AcSem.bestvalues{1}(:,AA)));
            %extremeBestDeviance(AA,2)=max(max(Deviance.Acoustic.bestvalues{1}(:,AA)),max(Deviance.AcSem.bestvalues{1}(:,AA)));
            legendcolor{AA,1}=num2str(Alphas(AA));
            subplot(2,2,4)
            hold on
            H(AA,3)=plot(1-R2.SS_Acoustic(:,AA)./R2.SS_T(:,AA),1-R2.SS_AcSem(:,AA)./R2.SS_T(:,AA),sprintf('%s.',Color(AA)), 'MarkerSize', 20);
        end
        subplot(2,2,1)
        legend(H(:,1),legendcolor{:,1})
        ylabel('AcSem Model R2')
        xlabel('AcSem Model Deviance')
        hold off
        subplot(2,2,2)
        legend(H(:,2),legendcolor{:,1})
        ylabel('Acoustic Model R2')
        xlabel('Acoustic Model Deviance')
        hold off
        subplot(2,2,4)
        legend(H(:,3),legendcolor{:,1})
        ylabel('AcSem Model R2')
        xlabel('Acoustic Model R2')
        axis([0 1 0 1])
        line([0 1],[0 1])
        hold off

        ss=subplot(2,2,3);
        MEAN_R2_AcSem = mean(1-R2.SS_AcSem./R2.SS_T,1);
        MEAN_R2_Acoustic = mean(1-R2.SS_Acoustic./R2.SS_T,1);
        MEAN_R2_Semantic = mean(1-R2.SS_Semantic./R2.SS_T,1);
        plot(1:NAlphas,MEAN_R2_AcSem,'r.','MarkerSize',20)
        hold on
        plot(MEAN_R2_Acoustic,'k.','MarkerSize',20)
        hold on
        plot(MEAN_R2_Semantic,'g.','MarkerSize',20)
        for AA=1:NAlphas
            hold on
            line([AA AA], [MEAN_R2_AcSem(AA)-std(1-R2.SS_AcSem(:,AA)./R2.SS_T(:,AA),1) MEAN_R2_AcSem(AA)+std(1-R2.SS_AcSem(:,AA)./R2.SS_T(:,AA),1)],'Color','r');
            hold on
            line([AA AA], [MEAN_R2_Acoustic(AA)-std(1-R2.SS_AcSem(:,AA)./R2.SS_T(:,AA),1) MEAN_R2_Acoustic(AA)+std(1-R2.SS_Acoustic(:,AA)./R2.SS_T(:,AA),1)],'Color','k');
            hold on
            line([AA AA], [MEAN_R2_Semantic(AA)-std(1-R2.SS_Semantic(:,AA)./R2.SS_T(:,AA),1) MEAN_R2_Semantic(AA)+std(1-R2.SS_Semantic(:,AA)./R2.SS_T(:,AA),1)],'Color','g');
        end
        hold on
        Pl=nan(3,1);
        Pl(1)=plot(1-sum(R2.SS_AcSem,1)./sum(R2.SS_T,1),'r*','MarkerSize',10);
        hold on
        Pl(2)=plot(1-sum(R2.SS_Acoustic,1)./sum(R2.SS_T,1),'k*','MarkerSize',10);
        hold on
        Pl(3)=plot(1-sum(R2.SS_Semantic,1)./sum(R2.SS_T,1),'g*','MarkerSize',10);
        set(ss,'XTick',[1:NAlphas])
        set(ss,'XTickLabel',Alphas)
        legend(Pl,'AcSem','Acoustic','Semantic','Location','SouthEast')
        ylabel('R-square')
        xlabel('Alphas')
        hold off
    end
    
    %% Effect of the exp output non-linearity?
    if NL_show
        for aa=1:length(Alphas)
            SpectroDim=size(Model.Acoustic.B{mm,aa});
            ModelAcoustic_local=repmat(reshape(Model.Acoustic.B{mm,aa}, 1, SpectroDim(1)*SpectroDim(2))./x_std,size(Data.x_wholeset{mm},1),1);
            Predict_Ac_wholeset = sum(ModelAcoustic_local.*Data.x_wholeset{mm},2) + repmat(Model.Acoustic.B0{mm,aa},size(Data.x_wholeset{mm},1),1);
            figure(17)
            ss=subplot(2,3,aa);
            plot(Predict_Ac_wholeset,Data.y_wholeset{mm},'b*','MarkerSize',10)
            ylabel('Observed Spike count ')
            xlabel('Predicted Spike count wo output non-linearity')
            title(sprintf('Acoustic Model Alpha %f',Alphas(aa)))
            AXES=get(ss,'yLim');
            hold on
            plot(Predict_Ac_wholeset,exp(Predict_Ac_wholeset),'r.')
            set(ss,'yLim',AXES)
            legend('Data','exponential fit')

            ModelAcSem_local=repmat([Model.AcSem.Bsem{mm,aa}' reshape(Model.AcSem.Bspectro{mm,aa}, 1, SpectroDim(1)*SpectroDim(2))./x_std],size(Data.x_wholeset{mm},1),1);
            Predict_AcSem_wholeset = sum(ModelAcSem_local.*[SemBoost.*Data.X_voc_wholeset{mm} Data.x_wholeset{mm}],2) + repmat(Model.AcSem.B0{mm,aa},size(Data.x_wholeset{mm},1),1);
            figure(18)
            ss=subplot(2,3,aa);
            plot(Predict_AcSem_wholeset,Data.y_wholeset{mm},'b*','MarkerSize',10)
            ylabel('Observed Spike count wo output non-linearity')
            xlabel('Predicted Spike count')
            title(sprintf('AcSem Model Alpha %f',Alphas(aa)))
            AXES=get(ss,'yLim');
            hold on
            plot(Predict_AcSem_wholeset,exp(Predict_AcSem_wholeset),'r.')
            set(ss,'yLim',AXES)
            legend('Data','exponential fit')

            ModelSem_local=repmat(Model.Semantic.B{mm,aa}',size(Data.x_wholeset{mm},1),1);
            Predict_Sem_wholeset = sum(ModelSem_local.*Data.X_voc_wholeset{mm},2) + repmat(Model.Semantic.B0{mm,aa},size(Data.x_wholeset{mm},1),1);
            figure(19)
            ss=subplot(2,3,aa);
            plot(Predict_Sem_wholeset,Data.y_wholeset{mm},'b*','MarkerSize',10)
            ylabel('Observed Spike count wo output non-linearity')
            xlabel('Predicted Spike count')
            title(sprintf('Semantic Model Alpha %f',Alphas(aa)))
            AXES=get(ss,'yLim');
            hold on
            plot(Predict_Sem_wholeset,exp(Predict_Sem_wholeset),'r.')
            set(ss,'yLim',AXES)
            legend('Data','exponential fit')
        end
        pause()
    end
end

% According to this code, 53.3 would be a good value for lambda
10^mean(MinDevlog10Lambdas.Acoustic)






