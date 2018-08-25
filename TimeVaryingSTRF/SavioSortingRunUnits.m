OutDir = '/auto/tdrive/julie/k6/julie/matfile/ModMatSavio';
DoneDir = '/auto/tdrive/julie/k6/julie/matfile/ModMatSavio/DoneCells';
TrashDir = '/auto/tdrive/julie/k6/julie/matfile/ModMatSavio/TrashCells';
UCDir = '/auto/tdrive/julie/k6/julie/matfile/ModMatSavio/UCCells';
SlurmoutFiles = dir(fullfile(OutDir, 'slurm_out*'));
Fnum = length(SlurmoutFiles);
CellToRunAgain = cell(Fnum,3);
TR=0;
%SlurmoutFiles = dir(fullfile(OutDir, 'slurm_out*'));
for ff=1:Fnum
    fprintf('File %d/%d--------------------------------------------\n',ff, Fnum)
    [status, cmdout] = system(sprintf('tail %s', fullfile(OutDir, SlurmoutFiles(ff).name)));
    cmdout
    DataFilename = ['Models_GLMPoisson_' SlurmoutFiles(ff).name(20:end-15) '.mat'];
    try
        Res=load(DataFilename);
        
        
        % Check if all windows were calculated
        UnrunWindows = find(isnan(Res.Deviance.Acoustic.FitIndex));
        if isempty(UnrunWindows)
            fprintf('All windows calculated, placing file in Done folder\n')
            %all windows were calculated move the data and results to the done folder
            [status, ~] = system(sprintf('cp %s %s', fullfile(OutDir, SlurmoutFiles(ff).name),fullfile(DoneDir, SlurmoutFiles(ff).name)));
            if ~status
                [~] = system(sprintf('rm %s %s', fullfile(OutDir, SlurmoutFiles(ff).name)));
            end
            [status, ~] = system(sprintf('cp %s %s', fullfile(OutDir, DataFilename),fullfile(DoneDir, DataFilename)));
            if ~status
                [~] = system(sprintf('rm %s %s', fullfile(OutDir, DataFilename)));
            end
        else
            fprintf('%d windows not calculated, file enter ToDocell\n',length(UnrunWindows))
            
            % establish that list of cell to be done and windows to be
            % calculated
            TR = TR +1;
            CellToRunAgain{TR,1} = fullfile(OutDir, SlurmoutFiles(ff).name);
            CellToRunAgain{TR,2} = fullfile(OutDir, DataFilename);
            CellToRunAgain{TR,3} = UnrunWindows;
        end
    catch err
        if strcmp(err.identifier, 'MATLAB:load:cantReadFile')
            fprintf('File unreadable, file enter ToDocell\n')
            % establish that list of cell to be done and windows to be
            % calculated
            TR = TR +1;
            CellToRunAgain{TR,1} = fullfile(OutDir, SlurmoutFiles(ff).name);
            CellToRunAgain{TR,2} = fullfile(OutDir, DataFilename);
            CellToRunAgain{TR,3} = 14;
            [status, ~] = system(sprintf('cp %s %s', fullfile(OutDir, SlurmoutFiles(ff).name),fullfile(TrashDir, SlurmoutFiles(ff).name)));
            if ~status
                [~] = system(sprintf('rm %s %s', fullfile(OutDir, SlurmoutFiles(ff).name)));
            end
            [status, ~] = system(sprintf('cp %s %s', fullfile(OutDir, DataFilename),fullfile(TrashDir, DataFilename)));
            if ~status
                [~] = system(sprintf('rm %s %s', fullfile(OutDir, DataFilename)));
            end
        elseif strcmp(err.identifier, 'MATLAB:load:couldNotReadFile')
            fprintf('Aborted Cell, file enter ToDocell\n')
            TR = TR +1;
            CellToRunAgain{TR,1} = fullfile(OutDir, SlurmoutFiles(ff).name);
            CellToRunAgain{TR,2} = fullfile(OutDir, DataFilename);
            CellToRunAgain{TR,3} = 14;
            [status, ~] = system(sprintf('cp %s %s', fullfile(OutDir, SlurmoutFiles(ff).name),fullfile(TrashDir, SlurmoutFiles(ff).name)));
            if ~status
                [~] = system(sprintf('rm %s %s', fullfile(OutDir, SlurmoutFiles(ff).name)));
            end
            [status, ~] = system(sprintf('cp %s %s', fullfile(OutDir, DataFilename),fullfile(TrashDir, DataFilename)));
            if ~status
                [~] = system(sprintf('rm %s %s', fullfile(OutDir, DataFilename)));
            end
        else
            fprintf('Not Sure what''s going on, file enter To DO cell\n')
            TR = TR +1;
            CellToRunAgain{TR,1} = fullfile(OutDir, SlurmoutFiles(ff).name);
            CellToRunAgain{TR,2} = fullfile(OutDir, DataFilename);
            CellToRunAgain{TR,3} = 14;
            [status, ~] = system(sprintf('cp %s %s', fullfile(OutDir, SlurmoutFiles(ff).name),fullfile(TrashDir, SlurmoutFiles(ff).name)));
            if ~status
                [~] = system(sprintf('rm %s %s', fullfile(OutDir, SlurmoutFiles(ff).name)));
            end
            [status, ~] = system(sprintf('cp %s %s', fullfile(OutDir, DataFilename),fullfile(TrashDir, DataFilename)));
            if ~status
                [~] = system(sprintf('rm %s %s', fullfile(OutDir, DataFilename)));
            end
        end
    end
    
    %pause()
end

% % Time limit savio
% slurmstepd: *** JOB 475843 CANCELLED AT 2015-12-11T14:06:37 DUE TO TIME LIMIT on n0134.savio0 ***
% 
% % Problem of file transfer but data ok
% Error writing to output stream.
% 
% 
% 
% Error in SpectroSemantic_Neuro_model_glm_savio (line 182)
%     [status1]=system(CommandSystemtransfer1);
%     
% % Failure in the pool start
% 
%     Failed to start pool.
%         Error using parallel.internal.pool.InteractiveClient>iAddFilesAndPaths
%         (line 691)
%         Unable to read MAT-file
%         /global/home/users/jelie/.matlab/local_cluster_jobs/R2015a/Job578.in.mat.
%         File might be corrupt.
%         
%         
%   %Same thing but different message
%   667)
%     Failed to initialize the interactive session.
%         Error using
%         parallel.internal.pool.InteractiveClient>iThrowIfBadParallelJobStatus
%         (line 768)
%         The interactive communicating job failed with no message.
%         
%         
%    %or:
%         
%    Error using parallel.internal.pool.InteractiveClient>iThrowWithCause (line
%     667)
%     Failed to start pool.
%         Error using parallel.Cluster/createCommunicatingJob (line 92)
%         The storage metadata file is corrupt. Please delete all files in the
%         JobStorageLocation and try again.



%% Renaming files
OutDir = '/auto/tdrive/julie/k6/julie/matfile/ModMatSavio';
cd /auto/tdrive/julie/k6/julie/matfile/ModMatSavio
ModelFiles = dir(fullfile(OutDir, 'Models_GLM*'));
for ff=1:length(ModelFiles)
    load(ModelFiles(ff).name);
    [pa,fi]=fileparts(MatfilePath);
    Newname = ['Models_GLMPoisson_' fi(10:end) '.mat'];
    [status, cmdout]=system(['mv ' ModelFiles(ff).name ' ' Newname])
end
    
    