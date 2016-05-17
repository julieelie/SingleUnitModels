function [Mat]=plotCellPerfAsMat(Signif, LocalCells, NumMod)
Mat = nan(length(LocalCells), NumMod,3);
for mm=1:NumMod
    % Id cells for which Ceil=Floor (white)
    LC=find(~Signif.CeilvsFl(LocalCells,mm));
    for ii=1:length(LC)
        Mat(LC(ii),mm,:)=[1 1 1];
    end
    
    % Id Cells for which Ceil>Floor and none of the two are better than
    % floor (Grey)
    LC=intersect(find(Signif.CeilvsFl(LocalCells,mm)),intersect(find(~Signif.AcvsFl(LocalCells,mm)),find(~Signif.SemvsFl(LocalCells,mm))));
    for ii=1:length(LC)
        Mat(LC(ii),mm,:)=[0.8 0.8 0.8];
    end
   
    % Id Cells for which Ceil>Floor and both Sem and Acoustic models are better than
    % Floor (orange)
    LC=intersect(find(Signif.CeilvsFl(LocalCells,mm)),intersect(find(Signif.AcvsFl(LocalCells,mm)),find(Signif.SemvsFl(LocalCells,mm))));
    for ii=1:length(LC)
        Mat(LC(ii),mm,:)=[1 0.5 0];
    end
    
    
    % Id Cells for which Ceil>Floor and only Ac is better than floor (yellow)
    LC=intersect(find(Signif.CeilvsFl(LocalCells,mm)),intersect(find(Signif.AcvsFl(LocalCells,mm)),find(~Signif.SemvsFl(LocalCells,mm))));
    for ii=1:length(LC)
        Mat(LC(ii),mm,:)=[1 1 0];
    end
    
    % Id Cells for which Ceil>Floor and only Sem is better than floor (red)
    LC=intersect(find(Signif.CeilvsFl(LocalCells,mm)),intersect(find(~Signif.AcvsFl(LocalCells,mm)),find(Signif.SemvsFl(LocalCells,mm))));
    for ii=1:length(LC)
        Mat(LC(ii),mm,:)=[1 0 0];
    end
    
end
end
