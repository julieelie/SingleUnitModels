OutputDir_final=fullfile('/auto','tdrive','julie','k6','julie','matfile','ModMatSavio');
OutputDir_local='/global/home/users/jelie/MatFiles/ModMat';
OutputDir_localscratch='/global/scratch/jelie/MatFiles/ModMat';
cd /global/scratch/jelie/MatFiles/ModMat
SavioFiles = dir(fullfile(OutputDir_localscratch, 'Models_GLMPoisson_*.mat'));
NF = length(SavioFiles);
for ff=1:NF
    fprintf('file %d/%d\n',ff,NF)
    [status, cmdout] = system(['ls -lt -h ' SavioFiles(ff).name])
    pause()
    StartIndex = regexp(cmdout, 'ucb');
    StartIndex = StartIndex(1)+3;
    EndIndex = regexp(cmdout, ' ');
    EndIndex = EndIndex(EndIndex>StartIndex);
    EndIndex = EndIndex(1)-1;
    if str2num(cmdout(StartIndex:EndIndex))==0
        fprintf('******Local File is empty, remove\n')
        system(['rm ' SavioFiles(ff).name],'-echo');
        system(['rm slurm_out*' SavioFiles(ff).name(19:end-4) '*.txt'],'-echo');
    else
        try
            LocalFile = load(SavioFiles(ff).name)
        catch ME
            fprintf('****error Loading the local file Removing the Local File\n')
            system(['rm ' SavioFiles(ff).name],'-echo');
            system(['rm slurm_out*' SavioFiles(ff).name(19:end-4) '*.txt'],'-echo');
            clear LocalFile
        end
        
        
        if ~exist('LocalFile', 'var')
            fprintf('*******Problem loading local file, likely corrupted, remove\n')
            system(['rm ' SavioFiles(ff).name],'-echo');
            system(['rm slurm_out*' SavioFiles(ff).name(19:end-4) '*.txt'],'-echo');
        else
            if ~(isfield(LocalFile, 'Deviance') && isfield(LocalFile, 'LL') && isfield(LocalFile, 'LambdaChoice') && isfield(LocalFile, 'Model') && isfield(LocalFile, 'PropVal') && isfield(LocalFile, 'Data') && isfield(LocalFile, 'Wins'))
                fprintf('******Local File is not complete, remove\n')
                system(['rm ' SavioFiles(ff).name],'-echo');
                system(['rm slurm_out*' SavioFiles(ff).name(19:end-4) '*.txt'],'-echo');
                clear LocalFile
            else
                fprintf('LocalFile has all the required fields of Data\n')
                Filecmd2 = ['ssh dtn.brc scp TheunissenLab:' OutputDir_final '/' SavioFiles(ff).name(1:end-4) '.mat' ' ' OutputDir_localscratch '/' SavioFiles(ff).name(1:end-4) '_old.mat'];
                [status,~] = system(Filecmd2, '-echo');
                if ~status
                    fprintf('Found some old data for this unit\n')
                    DoneCalc = load([OutputDir_localscratch '/' SavioFiles(ff).name(1:end-4) '_old.mat']);
                    if isfield(DoneCalc, 'Deviance') && isfield(DoneCalc, 'LL') && isfield(DoneCalc, 'LambdaChoice') && isfield(DoneCalc, 'Model') && isfield(DoneCalc, 'PropVal') && isfield(DoneCalc, 'Data') && isfield(DoneCalc, 'Wins')
                        fprintf('Old Data seems ok\n')
                        % Determine where the calculus stopped for old file
                        UnrunWindows_Old =[];
                        modNum = length(DoneCalc.Model.MeanSpectroStim);
                        for ww=1:modNum
                            if isempty(DoneCalc.Model.MeanSpectroStim{ww})
                                UnrunWindows_Old = [UnrunWindows_Old ww];
                            end
                        end
                        StartWin_Old = min(UnrunWindows_Old);
                        % Determine where the calculus stopped for local file
                        UnrunWindows_Local =[];
                        for ww=1:modNum
                            if isempty(LocalFile.Model.MeanSpectroStim{ww})
                                UnrunWindows_Local = [UnrunWindows_Local ww];
                            end
                        end
                        StartWin_Local = min(UnrunWindows_Local);
                        if StartWin_Local>StartWin_Old
                            fprintf('********Keep that file %s it is better than the older one run %d compare to %d windows\n',SavioFiles(ff).name, StartWin_Local,StartWin_Old)
                            system(['rm ' OutputDir_localscratch '/' SavioFiles(ff).name(1:end-4) '_old.mat'],'-echo')
                            system(['mv ' SavioFiles(ff).name ' GoodMat2Keep/'],'-echo');
                            system(['mv slurm_out*' SavioFiles(ff).name(19:end-4) '*.txt' ' GoodMat2Keep/'],'-echo');
                        else
                            fprintf('*****remove that file, not better than the one on the cluster run only %d compare to %d windows on the cluster\n', StartWin_Local,StartWin_Old)
                            system(['rm ' SavioFiles(ff).name],'-echo');
                            system(['rm slurm_out*' SavioFiles(ff).name(19:end-4) '*.txt'],'-echo');
                            system(['rm ' OutputDir_localscratch '/' SavioFiles(ff).name(1:end-4) '_old.mat'],'-echo')
                        end
                    else
                       fprintf('Old Data does not have all the required fields of Data\n')
                       fprintf('********Keep the new file %s it is better than the older one\n',SavioFiles(ff).name)
                       system(['rm ' OutputDir_localscratch '/' SavioFiles(ff).name(1:end-4) '_old.mat'],'-echo')
                       system(['mv ' SavioFiles(ff).name ' GoodMat2Keep/'],'-echo');
                       system(['mv slurm_out*' SavioFiles(ff).name(19:end-4) '*.txt' ' GoodMat2Keep/'],'-echo');
                    end
                    clear LocalFile
                else
                    fprintf('*******Could not find some old data for this file, keep the new file %s\n',SavioFiles(ff).name)
                            system(['mv ' SavioFiles(ff).name ' GoodMat2Keep/'],'-echo');
                            system(['mv slurm_out*' SavioFiles(ff).name(19:end-4) '*.txt' ' GoodMat2Keep/'],'-echo');
                    clear LocalFile
                end
            end
        end
        
    end
    clear LocalFile
end
            
            
            