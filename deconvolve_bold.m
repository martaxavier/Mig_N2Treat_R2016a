% This script estimates the BOLD events and the HRF
% parameters of each subject's HRF and deconvolves
% their BOLD signal accordingly. 

clear all
close all

set(groot, 'defaultFigureUnits','normalized');
set(groot, 'defaultFigurePosition',[0 0 1 1]);
set(0,'DefaultFigureVisible','off');

path='/home/mxavier/Mig_N2Treat/MATLAB'; cd(path)

% ------------------------------------------------------------
% Specify input and output data 
% ------------------------------------------------------------ 

% Specify list of subjects  
subjects = string({'sub-patient002',...
    'sub-patient003','sub-patient006',...
    'sub-pilot011','sub-pilot015','sub-pilot018'});
n_subjects = length(subjects);

% Specify input options 
mode = 'voxel_wise'; % 'voxel_wise', 'time_series'
bold = 'avg_dmn';     % 'avg_dmn', 'ic_dmn' 
report = 0;

% Specify input and output data
data_bold_in = 'bold_preprocessed.txt';
data_bold_img_in = 'bold_img_preprocessed.txt';
data_bold_out = 'bold_preprocessed_deconv';
data_bold_img_out = 'bold_img_preprocessed_deconv';
if strcmp(mode,'time_series'); data_in = data_bold_in; 
data_out = data_bold_out; elseif strcmp(mode,'voxel_wise')
data_in = data_bold_img_in; data_out = data_bold_img_out; 
end

% Specify input and output directories
if strcmp(mode,'time_series')
    path_data_in = strcat('DERIVATIVES/',subjects,'/func/',bold,'/');
    path_data_out = strcat('DERIVATIVES/',subjects,'/func/',bold,'/');
    path_img_out = strcat('RESULTS/',subjects,'/func/',bold,'/');
elseif strcmp(mode,'voxel_wise')
    path_data_in = strcat('DERIVATIVES/',subjects,'/func/');
    path_data_out = strcat('DERIVATIVES/',subjects,'/func/');
    path_img_out = strcat('RESULTS/',subjects,'/func/');
end

n_xvoxs = 100;
n_yvoxs = 100;
n_zvoxs = 60;

% ------------------------------------------------------------
% Write report in case of report 
% ------------------------------------------------------------  

if report
    
    import mlreportgen.dom.*;
    import mlreportgen.report.*;

    % Create the report object 
    R = Report(strcat('RESULTS -', ...
        'BOLD DECONVOLUTION'),'pdf');
    R.Layout.Landscape = true;
    open(R)

    % Create report heading 
    H = Heading(1, 'TIME-SERIES');
    H.FontSize='32pt';
    H.Color='#00b050';
    H.FontFamilyName='Arial';
    add(R,H)
    
    for s = 1 : n_subjects
        
        H2 = Heading(2,subjects(s));
        H2.FontSize='24pt';
        H2.Color='#809ec2';
        H2.FontFamilyName='Arial';
        add(R,H2)
        
        % Add image of BOLD, BOLD
        % deconvolved, BOLD events
        source = fullfile(path_img_out(s),...
            strcat(data_bold_out,'.png'));
        I = FormalImage(source);
        I.ScaleToFit=true; add(R,I);
        
    end
   
    % Create report heading 
    H = Heading(1, 'VOXEL-WISE');
    H.FontSize='32pt';
    H.Color='#00b050';
    H.FontFamilyName='Arial';
    add(R,H)
    
    for s = 1 : n_subjects
        
        H2 = Heading(2,subjects(s));
        H2.FontSize='24pt';
        H2.Color='#809ec2';
        H2.FontFamilyName='Arial';
        add(R,H2)
        
        % Add histogram images
        H3 = Heading(3,'Histograms');
        H3.FontSize = '14pt';
        H3.Color = '#404040';
        H3.FontFamilyName = 'Arial';
        add(R,H3);
        
        % Height
        source = fullfile(path_img_out(s),...
            'histo_estimated_hrfs_height.png');
        I = FormalImage(source);
        I.ScaleToFit=true; add(R,I);
        
        % Time-to-peak
        source = fullfile(path_img_out(s),...
            'histo_estimated_hrfs_timetopeak.png');
        I = FormalImage(source);
        I.ScaleToFit=true; add(R,I);
        
        % FWHM
        source = fullfile(path_img_out(s),...
            'histo_estimated_hrfs_fwhm.png');
        I = FormalImage(source);
        I.ScaleToFit=true; add(R,I);
        
        % Add histogram images
        H3 = Heading(3,'Maps');
        H3.FontSize = '14pt';
        H3.Color = '#404040';
        H3.FontFamilyName = 'Arial';
        add(R,H3);
        
        % Height
        source = fullfile(path_img_out(s),...
            'map_estimated_hrfs_height.png');
        I = FormalImage(source);
        I.ScaleToFit=true; add(R,I);
        
        % Time-to-peak
        source = fullfile(path_img_out(s),...
            'map_estimated_hrfs_timetopeak.png');
        I = FormalImage(source);
        I.ScaleToFit=true; add(R,I);
        
        % FWHM
        source = fullfile(path_img_out(s),...
            'map_estimated_hrfs_fwhm.png');
        I = FormalImage(source);
        I.ScaleToFit=true; add(R,I);
        
    end
    
    close(R)
    
    return
    
end

% ------------------------------------------------------------
% Specify input parameters
% ------------------------------------------------------------  

% Input information 
fs = 4;
temporal_mask = [];

% Input basis function structure
xBF.name = 'Canonical HRF (time and dispersion derivatives)';  

xBF.len = 32;                            % length in seconds of basis
xBF.order = 1;                           % order of basis set

xBF.T = 16;                              % number of subdivisions of TR
                                         % microtime resolution 
xBF.TR = 1/fs;                           % scan repetition time (seconds)
xBF.T0 = fix(xBF.T/2);                   % first time bin 
xBF.dt = xBF.TR/xBF.T;                   % length of time bin (seconds)
                                         % xBF.dt = xBF.TR/xBF.T
                                         
xBF.AR_lag = 1;                          % AR(1) noise autocorrelation
xBF.thr = 1;                             % threshold for BOLD events 
 
min_onset_search = 2;                    % minimum delay allowed between
                                         % event and HRF onset (seconds)
max_onset_search = 10;                   % maximum delay allowed between 
                                         % event and HRF onset (seconds)

xBF.lag = fix(min_onset_search/xBF.dt)...% array of acceptable lags (bins)
    :fix(max_onset_search/xBF.dt);

% ------------------------------------------------------------
% Go through subjects 
% ------------------------------------------------------------  

best_delay = zeros(n_subjects,1);

for s = 1 : n_subjects
      
    if ~exist(char(path_data_out(s)), 'dir'); mkdir(char(path_data_out(s))); end
    if ~exist(char(path_img_out(s)), 'dir'); mkdir(char(path_img_out(s))); end

    data = dlmread(strcat(path_data_in(s),data_in));
    n_pnts = size(data,1);
    n_voxs = size(data,2);
    
    if strcmp(mode,'voxel_wise')
        dmn = dlmread(strcat(path_data_in(s),'dmn_mask.txt'));
        data_img = zeros([size(dmn) 3]);
    end
    
    % ------------------------------------------------------------
    % Estimate BOLD events and HRF 
    % ---------------------a---------------------------------------  

    %disp(strcat('Estimating HRF',...
    %' for', " ", subjects(s),' ...'));
    tic
    
    % Estimate the HRF basis function and BOLD events 
    % hrf_beta(1) - scale; hrf_beta(2) - onset; hrf_beta(3) - best lag 
    % bf - orthogonalized HRF basis function 
    % event_bold - estimated pseudo-events
    [hrf_beta, bf, event_bold] = ...
        rsHRF_estimation_temporal_basis(data,xBF,temporal_mask);

    if strcmp(mode,'time_series');n_events=length(event_bold{1,1});end

    % Scale HRF basis function 
    % by the scaling parameter
    hrfa = bf*hrf_beta(1:size(bf,2),:);

    % Cycle parfor here?
    % Estimate HRF parameters from the 
    % hrf_pars(1) - height
    % hrf_pars(2) - time to peak (derivative < 0.001)
    % hrf_pars(3) - width 
    
    hrf_pars = zeros(3,n_voxs);
    for v = 1 : n_voxs
        hrf_pars(:,v) = ...
            wgr_get_parameters(hrfa(:,v),xBF.TR/xBF.T);
    end 
    
    if strcmp(mode,'voxel_wise')
        data_img = squeeze(data_img);
        data_img(find(dmn),1) = hrf_pars(1,:);
        data_img(find(dmn),2) = hrf_pars(2,:);
        data_img(find(dmn),3) = hrf_pars(3,:);
    
        data_img = reshape(data_img,...
            [n_xvoxs n_yvoxs n_zvoxs 3]);
    end
    
    toc
    disp('Done');
    
    % ------------------------------------------------------------
    % Retrieve estimated lag 
    % ------------------------------------------------------------  

    % Estimated HRF lag (seconds)
    %hrf_lag = hrf_beta(end)*xBF.dt;
    hrf_lag = hrf_pars(2,:); 

    % Save best lag for current subject 
    if strcmp(mode,'time_series');best_delay(s) = hrf_lag;end

    % ------------------------------------------------------------
    % Perform deconvolution of the BOLD signal (default lag)
    % ------------------------------------------------------------  

    %disp(strcat('Performing deconvolution',...
    %    ' for', " ", subjects(s),' ...'));
    tic

    hrfa_TR = resample(hrfa,1,xBF.T);
    hrf=hrfa_TR;

    % Deconvolution using a Wiener (restoration) filter 
    H=fft(cat(1,hrf, zeros(n_pnts-size(hrf,1),n_voxs)));
    B=fft(data);
    data_deconv = ifft(conj(H).*B./(H.*conj(H)+.1*mean(H.*conj(H))));
    
    % 0 mean and standard deviation 1 
    data_deconv = zscore(data_deconv);
    
    % Save deconvolved BOLD sinal in output directory
    dlmwrite(char(strcat(path_data_out(s),data_out,'.txt')), ...
        zscore(data_deconv));
    
    toc
    disp('Done');
    
    deconv_file = strcat('deconv_',mode,'.mat');
    save(char(strcat(path_data_out(s),deconv_file)),...
        'hrf_beta','hrf_pars','event_bold');

    % ------------------------------------------------------------
    % Time-series: plot deconvolved data and save results
    % ------------------------------------------------------------  
    
    if strcmp(mode,'time_series')
        
        % Plot HRF
        figure; 
        time = (1:length(hrfa(:,1)))*xBF.TR/xBF.T;
        plot(time,hrfa(:,1),'Color','#0072BD');
        delay_msg = strcat(upper('delay:'), num2str(hrf_lag), ' s');
        text(time(end)-5,1,delay_msg,'FontSize', 12); hold on;
        xlabel('Time (s)'); title('Estimated HRF');
        saveas(gcf,fullfile(path_img_out(s),strcat('estimated_hrf.png')));

        event_plot=zeros(1,n_pnts);
        event_plot(event_bold{1,1})=1;
        event_plot(~event_plot)=nan;

        time = (1:n_pnts)*xBF.TR;

        % For the HRF with the default lag 
        figure;
        text(time(end)-50,-3.5,delay_msg,'FontSize', 12); hold on;
        plot(time,zscore(data),'LineWidth',0.7,'Color','#0072BD'); hold on;
        plot(time,zscore(data_deconv),'LineWidth',0.7,'Color','#77AC30');
        stem(time,event_plot,'Color','k','LineWidth',1);
        legend('BOLD','BOLD deconvolved','BOLD events','FontSize',14, ...
            'FontWeight','normal'); legend('boxoff');
        xlabel('Time (s)');
        saveas(gcf,fullfile(path_img_out(s),strcat(data_out,'.png')));
        
        % NOTE - zscore returns the z-score for each element of X 
        % such that columns of X are centered to have mean 0 and 
        % scaled to have standard deviation 1 
     
    % ------------------------------------------------------------
    % Voxel-wise: plot histograms of HRF pars within the DMN 
    % ------------------------------------------------------------  
    
    elseif strcmp(mode,'voxel_wise')
        
        zslice = 35;
        
        % Plot histogram of height
        fig = figure('Name','Height - histogram');
        fig.Position(3:4) = fig.Position(3:4)*5;
        histogram(hrf_pars(1,:));
        title('Height of the estimated HRFs of voxels in the DMN');
        ax = gca; ax.XGrid = 'off'; ax.YGrid = 'on'; 
        set(gca,'FontSize',16); 
        img_out = 'histo_estimated_hrfs_height';
        saveas(gcf,fullfile(char(path_img_out(s)),img_out),'png');  
        
        % Plot map of height 
        fig = figure('Name','Height - map');
        fig.Position(3:4) = fig.Position(3:4)*5;
        imagesc(squeeze(data_img(:,:,zslice,1))); xlabel('x'); ylabel('y')
        title(char(strcat('Height of the estimated HRFs of voxels', ...
            'in the DMN, axial slice (z =', num2str(zslice),')')));
        set(gca,'FontSize',16); colorbar; 
        axis equal; img_out = 'map_estimated_hrfs_height';
        saveas(gcf,fullfile(char(path_img_out(s)),img_out),'png');  
        
        % Plot histogram of time-to-peak
        fig = figure('Name','Time-to-peak'); 
        fig.Position(3:4) = fig.Position(3:4)*5;
        histogram(hrf_pars(2,:));
        title('Time-to-peak of the estimated HRFs of voxels in the DMN',...
        'FontWeight','Normal');
        ax = gca; ax.XGrid = 'off'; ax.YGrid = 'on'; 
        set(gca,'FontSize',16); 
        img_out = 'histo_estimated_hrfs_timetopeak';
        saveas(gcf,fullfile(char(path_img_out(s)),img_out),'png');       
 
         % Plot map of time-to-peak 
        fig = figure('Name','Time-to-peak - map');
        fig.Position(3:4) = fig.Position(3:4)*5;
        imagesc(squeeze(data_img(:,:,zslice,2))); xlabel('x'); ylabel('y')
        title(char(strcat('Time-to-peak of the estimated HRFs of voxels', ...
            'in the DMN, axial slice (z =', num2str(zslice),')')));
        set(gca,'FontSize',16); colorbar;
        axis equal; img_out = 'map_estimated_hrfs_timetopeak';
        saveas(gcf,fullfile(char(path_img_out(s)),img_out),'png');  
        
        % Plot histogram of FWHM
        fig = figure('Name','FWHM');
        fig.Position(3:4) = fig.Position(3:4)*5;
        histogram(hrf_pars(3,:));
        title('FWHM of the estimated HRFs of voxels in the DMN',...
        'FontWeight','Normal');
        ax = gca; ax.XGrid = 'off'; ax.YGrid = 'on'; 
        set(gca,'FontSize',16); 
        img_out = 'histo_estimated_hrfs_fwhm';
        saveas(gcf,fullfile(char(path_img_out(s)),img_out),'png');  
        
         % Plot map of FWHM 
        fig = figure('Name','FWHM - map');
        fig.Position(3:4) = fig.Position(3:4)*5;
        imagesc(squeeze(data_img(:,:,zslice,3))); xlabel('x'); ylabel('y')
        title(char(strcat('FWHM of the estimated HRFs of voxels', ...
            'in the DMN, axial slice (z =', num2str(zslice),')')));
        set(gca,'FontSize',16); colorbar;
        axis equal; img_out = 'map_estimated_hrfs_fwhm';
        saveas(gcf,fullfile(char(path_img_out(s)),img_out),'png');  
        
    end

end % finish looping through subjects

if strcmp(mode,'time_series')

    % Save table containing best_lag for each subject 
    subjects = subjects';
    deconv_delay = table(subjects,best_delay);
    save(strcat(path,'\deconv_delay.mat'),'deconv_delay');

end

