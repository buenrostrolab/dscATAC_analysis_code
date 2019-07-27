% Author: Jason Buenrostro, Harvard University
% Compare stim-response
clear; clc
addpath(genpath('/Users/jasonbuenrostro/Documents/MATLAB/bin'));

% new anlaysis to do:
rowLabels = readtextJDB('./sci_data/scScores.rowLabels.txt');
scores = dlmread('./sci_data/scScores.txt','\t',1,3);

% load TFs and coords
tCoords = dlmread('./sci_data/tsneControls.txt','\t',1,0);
tfScores = dlmread('./sci_data/chromVar_Zscores.noNA.txt','\t',1,1);
tfScores(find(tfScores==0)) = nan; % fix the na's

% check sort
tfNames = readtextJDB('./sci_data/chromVar_Zscores.colLabels.txt');tfNames = strsplitJDB('_',tfNames,3);
scores_names = readtextJDB('./sci_data/scScores.rowLabels.txt');
tfScores_names = readtextJDB('./sci_data/chromVar_Zscores.rowLabels.txt');tfScores_names = tfScores_names(2:end);
sum(strcmp(scores_names(:,1),tfScores_names(:)))

%% split
% get conditions
idxC = find(strfindJDB('Control',rowLabels(:,3)));
idxS = find(strfindJDB('LPS',rowLabels(:,3)));

% smooth ref
clear smo_ref
b1 = knnsearch(scores(idxC,:),scores(idxC,:),'Distance','correlation','K',20);
for i = 1:length(idxC);disp(i)
    smo_ref(i,:) = mean(scores(idxC(b1(i,:)),:));
end

% get top 20 most similar
b2 = knnsearch(scores(idxS,:),smo_ref,'Distance','correlation','K',20);

%%
% loop through and get diff
tfScores_norm = tfScores;topN = 20;
clear tf1 tf2
for i = 1:size(b1,1);disp(i)
    tf1(i,:) = tfScores_norm( idxC(i) ,:);
    tf2(i,:) = mean(tfScores_norm( idxS(b2(i,1:topN)) ,:));
end
diffTFs = tf2-tf1;

% plot diff
figure;plot(std(diffTFs),'.');text(1:length(tfNames),std(diffTFs),tfNames)
figure;plot(std(tf1),std(diffTFs),'.');ylabel('Variable after stim norm.');xlabel('Variable')
text(std(tf1),std(diffTFs),tfNames)

% plot TF
tf = 139;
figure;scatter(tCoords(sub,1),tCoords(sub,2),[],diffTFs(sub,tf),'.');title(tfNames(tf));
title([tfNames(tf) 'score after stim']);

%% smooth a somewhat noisy signal
% smoth for k nearest
clear diffTFs_smo
b3 = knnsearch(scores(idxC(1:size(diffTFs,1)),:),scores(idxC,:),'Distance','correlation','K',20);
for i = 1:size(diffTFs,1);disp(i)
    diffTFs_smo(i,:) = mean(diffTFs(b3(i,1:topN),:));
end

% plot diff
x = std(tf1);y = std(diffTFs_smo);
figure;plot(x,y,'.');ylabel('Variable after stim norm.');xlabel('Variable')
text(x(x>1.3 | y >.4),y(x>1.3 | y >.4),tfNames(x>1.3 | y >.4))
title('LPS')

% plot less noisy TF
tf = 139;
figure;scatter(tCoords(sub,1),tCoords(sub,2),[],diffTFs_smo(sub,tf),'.');title(tfNames(tf));
title([tfNames(tf) 'score after stim & smoothed']);

% save
%dlmwrite(['lps_diff_' date '.txt'],diffTFs,'\t')
%dlmwrite(['lps_diff_' date '.smooth.txt'],diffTFs_smo,'\t')
%diffTFs_smoS = dlmread('serum_diff_15-Oct-2018.smooth.txt');
