%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Example script for Dynamic Correlation Toolbox
%
% File created by Martin Lindquist 07/22/14
%
% Makes use of functions from the UCSD_Garch toolbox by Kevin Shepard (Please see license agreement)
%
% Before running this script, begin by adding the DC_toolbox and all its subdirectories to the Matlab path.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc
clear
resultdir='F:\3_PNAS_paranoia\simulation\simu_dcc_null';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create simulated data set 
% rng('default')
% Define data dimensions
custom_cm = cbrewer('seq','YlGnBu',10);
method='L1';
TR=1;
NoiseSD=1;
filename=[resultdir filesep 'nulldata_noisesd10'];
for s=1:10
p = 2;         % Number of nodes
T = 600;        % Numer of time points

% Generate null data
mu = zeros(p,1);
Sigma = [2 0; 0 3];
dat=mvnrnd(mu,Sigma,T);     
% Add a little gaussian noise
%  rng(10*(2*sub+1000))
maxdata=round(max(max(dat)));
mindata=round(min(min(dat)));
Noisedata=NoiseSD*randn(T,p);
Noisemax=max(max(Noisedata));
Noisemin=min(min(Noisedata));
kcoe=(maxdata-mindata)/(Noisemax-Noisemin);
NoisedataN=kcoe*(Noisedata);
dat = dat + NoisedataN*NoiseSD;
% Note the input data has dimensions T-by-p (time by #nodes)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fit DCC
% Slower, more accurate version
[Ct1 ] = DCC(dat);
subplot(5, 1, 1)
plot(0:T-1, squeeze(Ct1(1,2,:)), 'Color',custom_cm(s,:),'LineWidth',3);
plot(0:T-1, 0, 'Color',[0.9 0.1 0.1],'LineWidth',3);
grid on
ylabel('DCC')
ylim([-1,1]);
set(gca, 'YTick', [-1, -0.5, 0, 0.5,1], ...                             % Change the axes tick marks
        'YTickLabel', {'-1', '-0.5', '0', '0.5', '1'}, ...  %   and tick labels
        'XTick', [100 200 300 400 500 600], ...
        'XTickLabel', {'100','200','300','400','500','600'}, ...
        'TickLength', [0 0]);
set(gca, 'TickDir', 'in', 'Xgrid', 'on'); 
set(gca, 'FontName','Arial','FontSize',24,'LineWidth', 2);
set(gcf,'Position',[10 10 2500*0.9 1080*0.9]);
% title('DCC - dynamic correlation between nodes 1 and 2')
hold on
% Ct1 is the dynamic correlation matrix 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fit sliding-window correlations
wsize = 22;
[tmp_dFC]=pp_ReHo_dALFF_dFC_gift(dat,method,TR,wsize);
Ct2=zeros(2,2,T);
for wi=1:T-wsize
Ct2(:,:,wi+wsize-1)=sf_vec2mat(2,tmp_dFC(wi));
end
% Ct2 is the sliding window correlation matrix 
% Plot some of the results
subplot(5, 1, 2)
plot(0:T-1, squeeze(Ct2(2,1,:)), 'Color',custom_cm(s,:),'LineWidth',3);
grid on
ylabel('SWFC')
ylim([-1,1]);
set(gca, 'YTick', [-1, -0.5, 0, 0.5,1], ...                             % Change the axes tick marks
        'YTickLabel', {'-1', '-0.5', '0', '0.5', '1'}, ...  %   and tick labels
        'XTick', [100 200 300 400 500 600], ...
        'XTickLabel', {'100','200','300','400','500','600'}, ...
        'TickLength', [0 0]);
set(gca, 'TickDir', 'in', 'Xgrid', 'on'); 
set(gca, 'FontName','Arial','FontSize',24,'LineWidth', 2);
set(gcf,'Position',[10 10 2500*0.9 1080*0.9]);
% title('SWFC - dynamic correlation between nodes 1 and 2')
hold on

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fit MTD
MTDwsize = 4;
Ct3=coupling(dat,MTDwsize);
tmpp=squeeze(Ct3(2,1,:));
tmppN=2*(tmpp-min(tmpp))/(max(tmpp)-min(tmpp))-1;
% Ct3 is the MTD correlation matrix 
% Plot some of the results
subplot(5, 1, 3)
plot(0:T-1, tmppN, 'Color',custom_cm(s,:),'LineWidth',3);
grid on
ylabel('MTD')
ylim([-1,1]);
set(gca, 'YTick', [-1, -0.5, 0, 0.5,1], ...                             % Change the axes tick marks
        'YTickLabel', {'-1', '-0.5', '0', '0.5', '1'}, ...  %   and tick labels
        'XTick', [100 200 300 400 500 600], ...
        'XTickLabel', {'100','200','300','400','500','600'}, ...
        'TickLength', [0 0]);
set(gca, 'TickDir', 'in', 'Xgrid', 'on'); 
set(gca, 'FontName','Arial','FontSize',24,'LineWidth', 2);
set(gcf,'Position',[10 10 2500*0.9 1080*0.9]);
% title('MTD - dynamic correlation between nodes 1 and 2')
hold on

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fit FLS
mu = 100;
Ct4 = yuan_DynamicBC_fls_FC(dat,mu);
% tmpp=squeeze(Ct3(2,1,:));
% tmppN=2*(tmpp-min(tmpp))/(max(tmpp)-min(tmpp))-1;
% Ct3 is the MTD correlation matrix 
% Plot some of the results
subplot(5, 1, 4)
plot(0:T-1, squeeze(Ct4(2,1,:)), 'Color',custom_cm(s,:),'LineWidth',3);
grid on
ylabel('FLS')
ylim([-1,1]);
set(gca, 'YTick', [-1, -0.5, 0, 0.5,1], ...                             % Change the axes tick marks
        'YTickLabel', {'-1', '-0.5', '0', '0.5', '1'}, ...  %   and tick labels
        'XTick', [100 200 300 400 500 600], ...
        'XTickLabel', {'100','200','300','400','500','600'}, ...
        'TickLength', [0 0]);
set(gca, 'TickDir', 'in', 'Xgrid', 'on'); 
set(gca, 'FontName','Arial','FontSize',24,'LineWidth', 2);
set(gcf,'Position',[10 10 2500*0.9 1080*0.9]);
% title('FLS - dynamic correlation between nodes 1 and 2')
hold on
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fit KF
YKF(1,:,:)=dat';
pKF=6;
ucKF=0.03;
FKF = dynet_SSM_KF(YKF,pKF,ucKF);
for i=1:T
FKFR(:,:,i)=icatb_corrcov(squeeze(FKF.R(:,:,i)));
end
subplot(5, 1, 5)
plot(0:T-1, squeeze(FKFR(1,2,:)), 'Color',custom_cm(s,:),'LineWidth',3);
hold on
grid on
ylabel('GLKF')
ylim([-1,1]);
set(gca, 'YTick', [-1, -0.5, 0, 0.5,1], ...                             % Change the axes tick marks
        'YTickLabel', {'-1', '-0.5', '0', '0.5', '1'}, ...  %   and tick labels
        'XTick', [100 200 300 400 500 600], ...
        'XTickLabel', {'100','200','300','400','500','600'}, ...
        'TickLength', [0 0]);
set(gca, 'TickDir', 'in', 'Xgrid', 'on'); 
set(gca, 'FontName','Arial','FontSize',24,'LineWidth', 2);
set(gcf,'Position',[10 10 2500*0.9 1080*0.9]);
% title('FLS - dynamic correlation between nodes 1 and 2')
hold on
% FSTOK = dynet_SSM_STOK(YKF,pKF);
% for i=1:T
% FSTOKR(:,:,i)=icatb_corrcov(squeeze(FSTOK.R(:,:,i)));
% end
% subplot(4, 1, 2)
% plot(0:T-1, squeeze(FSTOKR(1,2,:)), 'Color',custom_cm(s,:),'LineWidth',3);
% hold on
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Fit HMM
% MTDwsize = 4;
% Ct3=[corrmat]=yuan_fun_hmm_main(data,T,TR,K,outdir);
% tmpp=squeeze(Ct3(2,1,:));
% tmppN=2*(tmpp-min(tmpp))/(max(tmpp)-min(tmpp))-1;
% % Ct3 is the MTD correlation matrix 
% % Plot some of the results
% subplot(4, 1, 3)
% plot(0:T-1, tmppN, 'Color',custom_cm(s,:),'LineWidth',3);
% grid on
% ylabel('MTD')
% ylim([-1,1]);
% set(gca, 'YTick', [-1, -0.5, 0, 0.5,1], ...                             % Change the axes tick marks
%         'YTickLabel', {'-1', '-0.5', '0', '0.5', '1'}, ...  %   and tick labels
%         'XTick', [50 100 150 200], ...
%         'XTickLabel', {'50','100','150','200'}, ...
%         'TickLength', [0 0]);
% set(gca, 'TickDir', 'in', 'Xgrid', 'on'); 
% set(gca, 'FontName','Arial','FontSize',24,'LineWidth', 2);
% set(gcf,'Position',[10 10 2500*0.9 1080*0.9]);
% % title('MTD - dynamic correlation between nodes 1 and 2')
% hold on
end

% filename=[resultdir filesep 'Nulldata_cov05'];
print(1,'-dtiffn','-r300',filename);
close(1)






