%CallBack function of the timer JobTimer
function [IndOrd] = jobchecker(obj,event,NJobs, MatNameToDo, MatfileToDo, JobParams,SlurmParams)
global JobTimer
SQ = evalc('system(''squeue -u julie'')');
NJobs_current = length(strfind(SQ,'julie'));
NJobs_ToStart = NJobs-NJobs_current;
if NJobs_ToStart
    % get current list of indices of files to run
    IndOrd = get(obj, 'UserData');
    for jj=1:NJobs_ToStart
        IndexFile = IndOrd(1);
        IndOrd(1)=[];
        JobParams.out = fullfile(SlurmParams.resultsDirectory,sprintf('slurm_out_%s_%%j.txt', MatNameToDo{IndexFile}));
        JobParams.err = JobParams.out;
        icmd = sprintf(SlurmParams.cmd, MatfileToDo{IndexFile});
        fprintf(1,'Calling slurm_sbatch with command %s\n',icmd);
        slurm_sbatch(icmd,JobParams);
        if isempty(IndOrd)
            fprintf(1,'DONE Starting all Jobs!!!\n');
            stop(JobTimer)
            break
        end
    end
    %Update changes to list of indices of files to run
    set(obj, 'UserData', IndOrd);
end
end