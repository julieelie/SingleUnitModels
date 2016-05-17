Res.Site = 'Site2_L1100R1450_e14_s0_ss1';
OutputDir_final=fullfile('/auto','tdrive','julie','k6','julie','matfile','ModMatSavio');
OutputDir_local='/global/home/users/jelie/MatFiles/ModMat';
calfilename_local=fullfile(OutputDir_local,'Models_GLMPoissonSite2_L1100R1450_e14_s0_ss1.mat');
calfilename_final=fullfile(OutputDir_final,'Models_GLMPoissonSite2_L1100R1450_e14_s0_ss1.mat');
outfilename_local=fullfile(OutputDir_local,['slurm_out*' Res.Site '*.txt']);
outfilename_final=fullfile(OutputDir_final,['slurm_out*' Res.Site '*.txt']);

Filecopied=0;
while Filecopied==0
    CommandSystemtransfer1=['scp ' calfilename_local ' TheunissenLab:' calfilename_final];
    CommandSystemtransfer2=['scp ' outfilename_local ' TheunissenLab:' OutputDir_final];
    fprintf(1,'Attempting transfer 1 with command: %s\n',CommandSystemtransfer1);
    [status1]=system(CommandSystemtransfer1);
    fprintf(1,'Attempting transfer 2 with command: %s\n',CommandSystemtransfer2);
    [status2]=system(CommandSystemtransfer2);
    if ~(status1 || status2)
        fprintf(1,'files correctly transfered');
        pwd
        system(['cp ' OutputDir_local '/JobToDoSavio/Ex*' Res.Site '*.txt ' OutputDir_local '/JobDoneSavio/'])
        fprintf(1,'moved Ex file');
        fprintf(1,'Now remove files');
        system(['rm ' calfilename_local])
        system(['rm ' OutputDir_local '/JobToDoSavio/Ex*' Res.Site '*.txt'])
        system(['rm ' outfilename_local])
        Filecopied=1;
    end
end
quit