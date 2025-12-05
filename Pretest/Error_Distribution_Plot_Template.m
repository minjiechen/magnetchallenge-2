close all; clear all; clc;


%% Set the style for the plots
set(0,'DefaultAxesFontName', 'LM Roman 12')
set(0,'DefaultAxesFontSize', 20)
set(0,'DefaultAxesFontWeight', 'Bold')
set(0,'DefaultTextFontname', 'LM Roman 12')
set(0,'DefaultTextFontSize', 20)
set(0,'DefaultTextFontWeight', 'Bold')

set(groot, 'defaultFigureColor', [1 1 1]); % White background for pictures
set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
set(groot, 'defaultTextInterpreter', 'latex');
set(groot, 'defaultLegendInterpreter', 'latex');
% It is not possible to set property defaults that apply to the colorbar label. 
% set(groot, 'defaultAxesFontSize', 12);
% set(groot, 'defaultColorbarFontSize', 12);
% set(groot, 'defaultLegendFontSize', 12);
% set(groot, 'defaultTextFontSize', 12);
set(groot, 'defaultAxesXGrid', 'on');
set(groot, 'defaultAxesYGrid', 'on');
set(groot, 'defaultAxesZGrid', 'on');
set(groot, 'defaultAxesXMinorTick', 'on');
set(groot, 'defaultAxesYMinorTick', 'on');
set(groot, 'defaultAxesZMinorTick', 'on');
set(groot, 'defaultAxesXMinorGrid', 'on', 'defaultAxesXMinorGridMode', 'manual');
set(groot, 'defaultAxesYMinorGrid', 'on', 'defaultAxesYMinorGridMode', 'manual');
set(groot, 'defaultAxesZMinorGrid', 'on', 'defaultAxesZMinorGridMode', 'manual');

red = [200 36 35]/255;
blue = [40 120 181]/255;
gray = [170 170 170]/255;

% Update Material

Material = 'N87';

% Load the predicted and measurement data.
pred = csvread(fullfile(Material, 'pred.csv'));
meas = csvread(fullfile(Material, 'meas.csv'));

datalength = 1000;          % Data segment total length.
Nseqs= length(pred(:,1));   % Number of sequences predicted.
ttl_segs = 3;               % Total number of segments that define the length of prediction. in this case, 10% 50% and 90%
ratios = [0.9, 0.5, 0.1];   % Set by the user anywhere between 0.01 to 0.99

%% Reading the complete measurement file

% Read the full measurement file.
file_name = fullfile(Material, [Material, '_Testing_True.h5']);
if isfile(file_name)
    B_true = h5read(file_name, '/B_seq', [1, 1], [inf,inf]);
    H_true = h5read(file_name, '/H_seq', [1, 1], [inf,inf]);
    Loss = h5read(file_name, '/Loss', [1, 1], [inf,inf]);
else
    disp('Hpy file not been generated yet, quit execution')
    return
end
H_true =  H_true';
B_true =  B_true';
Loss = Loss';

% Read the partial data file.
file_name = fullfile(Material, [Material, '_Testing_Padded.h5']);
if isfile(file_name)
    H_past = h5read(file_name, '/H_seq', [1, 1], [inf,inf]);
else
    disp('Hpy file not been generated yet, quit execution')
    return
end
H_past =  H_past';


%% Calculate the 95% of the RMSE for each various future sequence

for n = 1 : ttl_segs  % Loop through each prediction length data seperately

    B = B_true((n-1)*length(H_true)/ttl_segs + 1 : n*length(H_true)/ttl_segs,  length(H_true(1 , :))*(1 - ratios(n)) + 1 : end);
    prediction = pred((n-1)*length(H_true)/ttl_segs + 1 : n*length(H_true)/ttl_segs, 1 : length(H_true(1 , :))*ratios(n));
    measurement = meas((n-1)*length(H_true)/ttl_segs + 1 : n*length(H_true)/ttl_segs, 1 : length(H_true(1 , :))*ratios(n));
    Core_loss = Loss((n-1)*length(H_true)/ttl_segs + 1 : n*length(H_true)/ttl_segs);

    % relative rmse of the sequences 
    rms_H = rms((prediction - measurement),2);
    rmse_pm = rms_H ./ rms(measurement, 2) * 100;
    
    % Energy loss calculation
    dB = diff(B, 1, 2);   
    H_mid_p = (prediction(:, 1:end-1) + prediction(:, 2:end)) / 2;
    H_mid_m = (measurement(:, 1:end-1) + measurement(:, 2:end)) / 2;
    energy_loss_p = H_mid_p .* dB; 
    energy_loss_m = H_mid_m .* dB;
    
    % Total normalized energy error
    diff_ene =  energy_loss_p - energy_loss_m;
    ene_error = sum(diff_ene , 2) ./ abs(Core_loss) * 100;


    if (n == 1) || (n == 2) || (n == 3) % 10%, 50%, 90% known
        
        figure; %   Plot sequence relative error.
        
        histogram(abs(rmse_pm), 100, 'FaceColor', blue, 'FaceAlpha',1,'Normalization','probability'); hold on;
        y1 = 0.02;
        y = linspace(0,y1,50);
        plot(mean(prctile(abs(rmse_pm), 95))*ones(size(y)),y,'--','Color',red,'LineWidth',2);
        text(prctile(abs(rmse_pm),95),y1,['95-Prct=',num2str(prctile(rmse_pm,95),4),'\%'],'Color',red);

        plot(mean(abs(rmse_pm))*ones(size(y)),y,'--','Color',red,'LineWidth',2);
        text(mean(abs(rmse_pm)),y1+0.02,['Ave Error =',num2str(mean(abs(rmse_pm)),4),'\%'],'Color',red);

        grid on;
        set(gcf,'Position',[850,550,780,430])

        xlabel('Relative Error of the Sequence [\%]');
        ylabel('Ratio of Data Points'); 
        title(['Sequence Relative Error for ', Material])
        subtitle(['Avg=',num2str(mean(abs(rmse_pm)),4),'\%, ', ...
            '95-Prct=',num2str(prctile(abs(rmse_pm),95),3),'\%, ', ...
            'Max=',num2str(max(abs(rmse_pm)),4),'\%'], ...
            'FontSize',18)

        figure; %  Plot energy against total core loss
        
        histogram(abs(ene_error), 100, 'FaceColor', blue, 'FaceAlpha',1,'Normalization','probability'); hold on;
        y1 = 0.03;
        y = linspace(0,y1,50);
        plot(prctile(abs(ene_error), 95)*ones(size(y)),y,'--','Color',red,'LineWidth',2);
        text(prctile(abs(ene_error),95),y1,['95-Prct=',num2str(prctile(ene_error,95),4),'\%'],'Color',red);
        plot(mean(abs(ene_error))*ones(size(y)),y,'--','Color',red,'LineWidth',2);
        text(mean(abs(ene_error)),y1+0.02,['Ave Error =',num2str(mean(abs(ene_error)),4),'\%'],'Color',red);

        xlabel('Relative Energy Density [\%]');
        ylabel('Ratio of Data Points'); 

        grid on;
        set(gcf,'Position',[850,550,780,430])
        title(['Total Energy Normalized Relative Error for ', Material])
        subtitle(['Avg=',num2str(mean(abs(ene_error)),4),'\%, ', ...
            '95-Prct=',num2str(prctile(abs(ene_error),95),3),'\%, ', ...
            'Max=',num2str(max(abs(ene_error)),4),'\%'], ...
            'FontSize',18)
 
    end
end