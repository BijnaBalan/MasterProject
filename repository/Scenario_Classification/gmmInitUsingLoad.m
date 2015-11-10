clc;
clear all;
close all;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ProcessAll = 0;
NUM_OF_SCENARIOS = 3;
scearioArray = {'Corridor Following'; 'Dead End'; 'Obstacle Avoidance';'Scenario Unidentified'};
CORRIDOR_FOLLOWING = 1;
DEAD_END = 2;
OBSTACLE_AVOIDANCE = 3;
SCENARIO_UNIDENTFIED = 4;
numofMixtures = 3;
DATA_INDEX = 1;
COLOR_INDEX = 2;
ID_INDEX = 3;
MODEL_INDEX = 4;
OPTI_MODEL_INDEX = 5;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%if already model is present dont generate
if(~exist('GMMmodel.mat','file'))
    allDataSave = {};
    % extract the features
    load('CF.mat');
    allDataSave{1,DATA_INDEX} = feCombined';
    allDataSave{1,COLOR_INDEX} = [1 0 0];
    allDataSave{1,ID_INDEX} = 'Corridor Following';
    load('DE.mat');
    allDataSave{2,DATA_INDEX} = feCombined';
    allDataSave{2,COLOR_INDEX} = [0 1 0];
    allDataSave{2,ID_INDEX} = 'Dead End';
    load('OA.mat');
    allDataSave{3,DATA_INDEX} = feCombined';
    allDataSave{3,COLOR_INDEX} = [0 0 1];
    allDataSave{3,ID_INDEX} = 'Obstacle Avoidance';
    % training
    tdataFig = figure('Name','Training  Data');
    title('Comparison Feature values across Different Scenarios');
    xlabel('Entropy');
    ylabel('Edge Pixel Ratio');
    hold on;
    grid on;
    for irun = 1:NUM_OF_SCENARIOS
        data_mat = allDataSave{irun,DATA_INDEX};
        color = allDataSave{irun,COLOR_INDEX};
        [Priors, Mu, Sigma] = EM_init_kmeans(data_mat, numofMixtures);
        [Priors, Mu, Sigma] = EM(data_mat, Priors, Mu, Sigma);
        % model
        allDataSave{irun,MODEL_INDEX}= {Priors, Mu, Sigma};
        figp = plot(data_mat(1,:),data_mat(2,:));
        set(figp,'Color',color,'LineWidth',2);
    end
    legend('Corridor Following','Dead End','Obstacle Avoidance');
    hold off;
    grid off;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % plot GMM model for each scenarios
    fontsize=16;
    % idStr = allDataSave{irun,ID_INDEX};
    figure('Name','GMM - Across Different Scenarios');
    hold on;
    grid on;
    for irun = 1:NUM_OF_SCENARIOS
        data_mat = allDataSave{irun,DATA_INDEX};
        color = allDataSave{irun,COLOR_INDEX};
        model = allDataSave{irun,MODEL_INDEX};
        Mu = model{2};
        Sigma = model{3};
        plot(data_mat(1,:),data_mat(2,:),'Color',color);
        plotGMM(Mu, Sigma, color, 1);
        xlabel('Entropy','FontSize',fontsize);ylabel('Edge Pixel Ratio','FontSize',fontsize);
    end
    hold off;
    grid off;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % %optimised model
    nRep = 3;
    nbStates_min = 1;
    nbStates_max = 3;
    for irun=1:NUM_OF_SCENARIOS
        BIC_min = +inf;
        for j=nbStates_min:nbStates_max
            for k=1:nRep
                data_mat = allDataSave{irun,DATA_INDEX};
                [Priors, Mu, Sigma] = EM_init_kmeans(data_mat, j);
                [Priors, Mu, Sigma] = EM(data_mat, Priors, Mu, Sigma);
                
                Pxi = [];
                for m=1:j
                    Pxi(:,m) = gaussPDF(data_mat, Mu(:,m), Sigma(:,:,m));
                end
                
                Pxi = Pxi*Priors';
                LL = sum(log(Pxi));
                
                [d, n] = size(data_mat);
                nParam = d*(d+1)*j/2 + j-1 + j*d;
                
                BIC_tmp = -2 * (LL) + nParam * log(n);
                
                if(BIC_tmp < BIC_min)
                    BIC_min = BIC_tmp;
                    %optimized Model
                    allDataSave{irun,OPTI_MODEL_INDEX} = {Priors, Mu, Sigma};
                end
            end
        end
    end
    save('GMMmodel.mat', 'allDataSave');
else
    % load the data
    load('GMMmodel.mat');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot modified GMM model for each scenarios
fontsize=16;
figure('Name','Modified GM Model used for identifying Different Scenarios');
hold on;
grid on;
for irun = 1:NUM_OF_SCENARIOS
    data_mat = allDataSave{irun,DATA_INDEX};
    color = allDataSave{irun,COLOR_INDEX};
    model = allDataSave{irun,OPTI_MODEL_INDEX};
    Mu = model{2};
    Sigma = model{3};
    plot(data_mat(1,:),data_mat(2,:),'Color',color);
    plotGMM(Mu, Sigma, color, 1);
    xlabel('Entropy','FontSize',fontsize);ylabel('Edge Pixel Ratio','FontSize',fontsize);
end
hold off;
grid off;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Test against the obtained model
% now i have 3 models now
% now load the test image data set and perform the GMM.
testImageLocStruct(1).Location = 'G:\project\evorob\groups\vision\workingDirectories\Bijna\Ind_task\Test_Data\CF\New\set1';
testImageLocStruct(2).Location = 'G:\project\evorob\groups\vision\workingDirectories\Bijna\Ind_task\Test_Data\DE\New\set1';
testImageLocStruct(3).Location = 'G:\project\evorob\groups\vision\workingDirectories\Bijna\Ind_task\Test_Data\CF\New\set2';
testImageLocStruct(4).Location = 'G:\project\evorob\groups\vision\workingDirectories\Bijna\Ind_task\Test_Data\OA\New\set1';
testImageLocStruct(5).Location = 'G:\project\evorob\groups\vision\workingDirectories\Bijna\Ind_task\Test_Data\DE\New\set2';
testImageLocStruct(6).Location = 'G:\project\evorob\groups\vision\workingDirectories\Bijna\Ind_task\Test_Data\CF\New\set3';
testImageLocStruct(7).Location = 'G:\project\evorob\groups\vision\workingDirectories\Bijna\Ind_task\Test_Data\OA\New\set2';
% testImageLocStruct(8).Location = 'G:\project\evorob\groups\vision\workingDirectories\Bijna\Ind_task\Test_Data\DE\set3';
% testImageLocStruct(9).Location = 'G:\project\evorob\groups\vision\workingDirectories\Bijna\Ind_task\Test_Data\CF\set4';
% testImageLocStruct(10).Location = 'G:\project\evorob\groups\vision\workingDirectories\Bijna\Ind_task\Test_Data\OA\set4';
% testImageLocStruct(11).Location = 'G:\project\evorob\groups\vision\workingDirectories\Bijna\Ind_task\Test_Data\DE\set4';
% testImageLocStruct(12).Location = 'G:\project\evorob\groups\vision\workingDirectories\Bijna\Ind_task\Test_Data\CF\set6';
% testImageLocStruct(13).Location = 'G:\project\evorob\groups\vision\workingDirectories\Bijna\Ind_task\Test_Data\OA\set6';
% testImageLocStruct(14).Location = 'G:\project\evorob\groups\vision\workingDirectories\Bijna\Ind_task\Test_Data\DE\set5';
numTestFolders = size(testImageLocStruct,2);
for irun =1:numTestFolders
    scenarioDetected = SCENARIO_UNIDENTFIED;
    testImageLoc = testImageLocStruct(irun).Location;
    featurestest = (extractImFeat4SceneClassification(testImageLoc,ProcessAll))';
    numofImages = size(featurestest,2);
    LL_max = -inf;
    for i=1:NUM_OF_SCENARIOS
        modelAns =  allDataSave{i,OPTI_MODEL_INDEX};
        Priors = cell2mat(modelAns(1,1));
        Mu = cell2mat(modelAns(1,2));
        Sigma = cell2mat(modelAns(1,3));
        Pxi = [];
        for m=1:size(Mu,2)
            Pxi(:,m) = gaussPDF(featurestest, Mu(:,m), Sigma(:,:,m));
        end
        Pxi = Pxi*Priors';
        LL = sum(log(Pxi));
        if(LL > LL_max)
            LL_max = LL;
            scenarioDetected = i;
        end
    end
    disp('=================================================================================');
    disp([num2str(irun) '. Scenario Detected - '  cell2mat(scearioArray(scenarioDetected))]);
    disp([num2str(irun) '. Number of Images - '  num2str(numofImages)]);
    disp('=================================================================================');
end

