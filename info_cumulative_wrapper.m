function [Model] = info_cumulative_wrapper(ParamModel,SWITCH,Model,firstwin,mm,x_stim_indices_wholeset, Stim_local,FolderTempInfStorage)
fprintf('Pre-process data for cumulative information calculations\n')
X_stim_indices_wholeset = x_stim_indices_wholeset(firstwin:mm);
            mm_local=0;
            %Pre-process adat for a parfor loop
            P_YgivenS_allModel = cell(10,1);
            if ParamModel.ModelChoice(1) && ~SWITCH.AllAlpha
                % ACoustic Model
                mm_local=mm_local+1;
                P_YgivenS_allModel{mm_local}=Model.Acoustic.P_YgivenS_all1(firstwin:mm,1);
            end
            if ParamModel.ModelChoice(2) && ~SWITCH.AllAlpha
                % Semantic Model
                mm_local=mm_local+1;
                P_YgivenS_allModel{mm_local}=Model.Semantic.P_YgivenS_all1(firstwin:mm,1);
            end
            if ParamModel.ModelChoice(4) && ~SWITCH.AllAlpha
                % AcSemAc
                mm_local=mm_local+1;
                P_YgivenS_allModel{mm_local}=Model.AcSemAc.P_YgivenS_all1(firstwin:mm,1);
            end
            if ParamModel.ModelChoice(5) && ~SWITCH.AllAlpha
                % AcSemSem
                mm_local=mm_local+1;
                P_YgivenS_allModel{mm_local}=Model.AcSemSem.P_YgivenS_all1(firstwin:mm,1);
            end
            % Floor
            mm_local=mm_local+1;
            P_YgivenS_allModel{mm_local}=Model.Floor.P_YgivenS_all1(firstwin:mm,1);
            %Ceiling
            mm_local=mm_local+1;
            P_YgivenS_allModel{mm_local}=Model.Ceiling.P_YgivenS_all1(firstwin:mm,1);

            % AR
            if SWITCH.AR
                mm_local=mm_local+1;
                P_YgivenS_allModel{mm_local}=Model.AR.P_YgivenS_all1(firstwin:mm,1);
            end

            Cum_info = nan(mm_local,1);
            fprintf('Calculate Cumulative information for all models from win %d\n',firstwin)
            for modelrun=1:mm_local
                fprintf('Cumulative info %d/%d\n', modelrun, mm_local);
                Cum_info(modelrun) = info_cumulative_model_Calculus(P_YgivenS_allModel{modelrun},firstwin, mm,X_stim_indices_wholeset, Stim_local,modelrun,FolderTempInfStorage);
            end
               
            mm_local=0;
            %Post-process data for a parfor loop
            if ParamModel.ModelChoice(1) && ~SWITCH.AllAlpha
                fprintf('**CumInfo on Acoustic**\n')
                % ACoustic Model
                mm_local=mm_local+1;
                Model.Acoustic.(sprintf('cum_info%d',firstwin))(mm-firstwin+1,1)=Cum_info(mm_local);
            end
            if ParamModel.ModelChoice(2) && ~SWITCH.AllAlpha
                fprintf('**CumInfo on Semantic**\n')
                % Semantic Model
                mm_local=mm_local+1;
                Model.Semantic.(sprintf('cum_info%d',firstwin))(mm-firstwin+1,1)=Cum_info(mm_local);
            end
            if ParamModel.ModelChoice(4) && ~SWITCH.AllAlpha
                fprintf('**CumInfo on AcSemAc**\n')
                % AcSemAc
                mm_local=mm_local+1;
                Model.AcSemAc.(sprintf('cum_info%d',firstwin))(mm-firstwin+1,1)=Cum_info(mm_local);
            end
            if ParamModel.ModelChoice(5) && ~SWITCH.AllAlpha
                fprintf('**CumInfo on AcSemSem**\n')
                % AcSemSem
                mm_local=mm_local+1;
                Model.AcSemSem.(sprintf('cum_info%d',firstwin))(mm-firstwin+1,1)=Cum_info(mm_local);
            end
            % Floor
            fprintf('**CumInfo on Floor**\n')
            mm_local=mm_local+1;
            Model.Floor.(sprintf('cum_info%d',firstwin))(mm-firstwin+1,1)=Cum_info(mm_local);
            %Ceiling
            fprintf('**CumInfo on Ceiling**\n')
            mm_local=mm_local+1;
            Model.Ceiling.(sprintf('cum_info%d',firstwin))(mm-firstwin+1,1)=Cum_info(mm_local);
            % AR
            if SWITCH.AR
                fprintf('**CumInfo on AR**\n')
                mm_local=mm_local+1;
                Model.AR.(sprintf('cum_info%d',firstwin))(mm-firstwin+1,1)=Cum_info(mm_local);
            end