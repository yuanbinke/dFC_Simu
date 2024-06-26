
clc
clear
resultdir='F:\3_PNAS_paranoia\simulation\simu_dcc_periodic';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create simulated data set 
% rng('default')
% Define data dimensions
custom_cm = cbrewer('seq','YlGnBu',10);
method='L1';
TR=1;
Fs=TR;
NoiseSD=0.3;
% for s=5:10
p = 2;         % Number of nodes
T = 600;        % Numer of time points
% Generate null data
mu = zeros(p,1);
Sigma = zeros(2,2);
Sigma(1,1)=2;
Sigma(2,2)=3;
Sigma(2,1)=0;
Sigma(1,2)=0;

dat=mvnrnd(mu,Sigma,T);
maxdata=round(max(max(dat)));
mindata=round(min(min(dat)));
Noisedata=NoiseSD*randn(T,p);
Noisemax=max(max(Noisedata));
Noisemin=min(min(Noisedata));
kcoe=(maxdata-mindata)/(Noisemax-Noisemin);
NoisedataN=kcoe*(Noisedata);
dat2=dat + NoisedataN*NoiseSD;

paddedLength = 2^nextpow2(T); %2^nextpow2(sampleLength);
    theFreqSeries=fft(dat2,paddedLength);
    theSampleFreq=1/TR ;
    theFreqPrecision=theSampleFreq/paddedLength;
    theFreqLim =[theFreqPrecision: theFreqPrecision :theSampleFreq/2];
    theXLim =[2,(paddedLength/2 +1)];	%don't Contain DC, because AFNI don't contain it in PowerSpectrum
    
    %Calcute the Power Spectrum
    theFreqSeries =abs(theFreqSeries([theXLim(1):theXLim(2)],:)); % Get the half's amplitude
    theFreqSeries(1:end,:) =theFreqSeries(1:end,:).^2 /T;%Don't containt the DC component because abs didn't make DC 2-times to its original amplitude , dawnsong 20070629
    %theFreqSeries(1) =theFreqSeries(1) /length(theTimeCourse);	% now process the DC component
    
    %Since we dropped half the FFT, we multiply mx by 2 to keep the same energy.
    % The DC component and Nyquist component, if it exists, are unique and should not
    % be mulitplied by 2.
    theFreqSeries(1:end-1,:) =theFreqSeries(1:end-1,:) *2;
    
subplot(3, 1, 1)
plot(0:T-1, dat2, 'LineWidth',3);
grid on
ylabel('Amplitude (a.u.)')
% ylim([-1,1]);
% set(gca, 'YTick', [-1, -0.5, 0, 0.5,1], ...                             % Change the axes tick marks
%         'YTickLabel', {'-1', '-0.5', '0', '0.5', '1'}, ...  %   and tick labels
%         'XTick', [100 200 300 400 500 600], ...
%         'XTickLabel', {'100','200','300','400','500','600'}, ...
%         'TickLength', [0 0]);
legend('TS1','TS2','Location','NorthEastOutside');
set(gca, 'TickDir', 'in', 'Xgrid', 'on'); 
set(gca, 'FontName','Arial','FontSize',24,'LineWidth', 2);
% title('DCC - dynamic correlation between nodes 1 and 2')
hold on

subplot(3, 1, 2)
plot(1:theXLim(2)-1, theFreqSeries(:,1),'LineWidth',3);
hold on
plot(1:theXLim(2)-1, theFreqSeries(:,2),'LineWidth',3);
hold on
grid on
ylabel('Power (f)')
xlabel('Frequency (f)')
legend('TS1','TS2','Location','NorthEastOutside');
xlim([0,515]);
% set(gca, 'YTick', [-1, -0.5, 0, 0.5,1], ...                             % Change the axes tick marks
%         'YTickLabel', {'-1', '-0.5', '0', '0.5', '1'}, ...  %   and tick labels
%         'XTick', [100 200 300 400 500 600], ...
%         'XTickLabel', {'100','200','300','400','500','600'}, ...
%         'TickLength', [0 0]);
set(gca, 'XTick', [21 103 205 307 410 512], ...
        'XTickLabel', {'0.01','0.1','0.2','0.3','0.4','0.5'}, ...
        'TickLength', [0 0]);
set(gca, 'TickDir', 'in', 'Xgrid', 'on'); 
set(gca, 'FontName','Arial','FontSize',24,'LineWidth', 2);
% title('DCC - dynamic correlation between nodes 1 and 2')
hold on

subplot(3, 1, 3)
plot(0:T-1, zeros(600,1), 'LineWidth',3);
grid on
ylabel('Coefficient')
hold on
MTDwsize = 5;
Ct3=coupling(dat2,MTDwsize);
tmpp1=squeeze(Ct3(1,2,:));
% tmppN1=2*(tmpp1-min(tmpp1))/(max(tmpp1)-min(tmpp1))-1;
% Ct3 is the MTD correlation matrix 
% Plot some of the results
% RealTMTD(:,s)=tmppN1;
plot(0:T-1, tmpp1, 'LineWidth',3);
xlim([0,600]);
grid on
legend('Truth','MTD','Location','NorthEastOutside');
% ylabel('Coefficient')
% ylim([-1,1]);
% set(gca, 'YTick', [-1, -0.5, 0, 0.5,1], ...                             % Change the axes tick marks
%         'YTickLabel', {'-1', '-0.5', '0', '0.5', '1'}, ...  %   and tick labels
%         'XTick', [100 200 300 400 500 600], ...
%         'XTickLabel', {'100','200','300','400','500','600'}, ...
%         'TickLength', [0 0]);
set(gca, 'TickDir', 'in', 'Xgrid', 'on'); 
set(gca, 'FontName','Arial','FontSize',24,'LineWidth', 2);
set(gcf,'Position',[10 10 2500*0.8 1080*0.8]);
% title('DCC - dynamic correlation between nodes 1 and 2')
hold on


cd(resultdir)
save demofft_null.mat
filename=[resultdir filesep 'null_fft'];
print(1,'-dtiffn','-r300',filename);
close 1