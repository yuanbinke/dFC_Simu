clear all; close all; clc

Simsubs=100;
% SimTRs=10000;
SimROIs=10;
NoiseSD=0.6;

Rootdesdir='I:\ISDCC\Newsimtb_0918';
Desdir=[Rootdesdir filesep 'dFC_nozscore_1TR_simu' filesep 'NoiseSD_FIXHRF_pu01_au05_Poisson2345_NoiseSD' num2str(NoiseSD) filesep  'Simu10' filesep 'tc'];mkdir(Desdir)
% seed = 3571; % choose seed for repeatable simulation
% simtb_rand_seed(seed);

rng('default')
rng(1000) % to enable repeated generation of the same simulation

% number of components
nC = SimROIs;

% TR
TR = 2;

% number of different connectivity states
nStates = 4;

% probability of unique events
pU = 0.1;

% amplitude of unique events (relative to module-specific events)
aU = 0.5;

% probability of state specific events 
pState = .5;

%Module membership for each state
ModMem = zeros(nC,nStates);

% Number of event types (i.e., number of different modules)
nE = 3;

% Modules are formed by sharing of events.  In each state there are up to
% nE different modules. The matrix ModMem stores the membership of each
% component to a module.Note that module 2 in state 2 has nothing to do with 
% module 2 in the other states, it's just an index.
% Negative numbers indicate activity is negatively related to the events in
% the given module.
ModMem(1,:) = [2   -2   3    2];
ModMem(2,:) = [2   -2   3    2];
ModMem(3,:) = [2   -2   3    2];
ModMem(4,:) = [-2  3   3    2];
ModMem(5,:) = [-2  3   2    2];
ModMem(6,:) = [-2  3   2    2];
ModMem(7,:) = [-2  2   2    1];
ModMem(8,:) = [1   2    1   1];
ModMem(9,:) = [1   2    1   -2];
ModMem(10,:)= [1   2    1   -2];

% The sojourn distributions in data sets one and
%two were Poisson distributed with λ ¼ 6; 12; 20; and 30 for states 1,2,3
%and 4 respectively
SdwellS1=poissrnd(20,1,40);
% pdf_values = poisspdf(SdwellS1, 6);
SdwellS2=poissrnd(30,1,40);
SdwellS3=poissrnd(40,1,40);
SdwellS4=poissrnd(50,1,40);

%plot
histfitb(SdwellS1)
hold on
histfity(SdwellS2)
hold on
histfitg(SdwellS3)
hold on
histfitr(SdwellS4)
hold on
ylabel('Sojourn density')
xlabel('Time in state (s)')
set(gca, 'XTick', [10 20 30 40 50 60 70], ...                             % Change the axes tick marks
         'XTickLabel', {'10', '20', '30','40','50','60','70'}, ...  %   and tick labels
         'TickLength', [0 0]); grid on
ylim([0 0.12])
xlim([0 70])
set(gca, 'FontName','Arial','FontSize',24,'LineWidth', 3.5);
set(gcf,'Position',[10 10 2560*0.3 1440*0.3]);
box on
legend([1 3 5 7], {'S1','S2','S3','S4'},'FontSize',26);
cd(Rootdesdir)
filename=[Rootdesdir filesep 'SimuDensity_poisson_20304050_NoiseSD062'];
% print(1,'-dtiffn','-r300',filename);
% close 1
% SdwellS5=poissrnd(40,1,100);
% SdwellS6=poissrnd(50,1,100);
Sdwell1=[SdwellS1;ones(1,40)]';
Sdwell2=[SdwellS2;ones(1,40)*2]';
Sdwell3=[SdwellS3;ones(1,40)*3]';
Sdwell4=[SdwellS4;ones(1,40)*4]';
% Sdwell5=[SdwellS5;ones(1,100)*5]';
% Sdwell6=[SdwellS6;ones(1,100)*6]';
SdwellRaw=[Sdwell1;Sdwell2;Sdwell3;Sdwell4];
% state vector
AAA=randperm(length(SdwellRaw));
Sdwell=SdwellRaw(AAA,1);
Sorder=SdwellRaw(AAA,2);

% number of time points
nT = sum(Sdwell);
SimsubTCs=zeros(nT,Simsubs,SimROIs);
%% Create figure of the connectivity matrix for each state
% F = figure('color','w','Name', 'sim_neural_connectivity');
CMAll=zeros(nC,nC,4);
for ii = 1:nStates
%     subplot(1,nStates,ii)
    CM = zeros(nC,nC);
    for jj = 1:nC
        for kk = 1:nC
            if ModMem(jj,ii) == ModMem(kk,ii)
                CM(jj,kk) = 1;
            elseif abs(ModMem(jj,ii)) == abs(ModMem(kk,ii))
                CM(jj,kk) = -1;
            else
                CM(jj,kk) = 0;
            end
        end
    end
    CMAll(:,:,ii)=CM;
%     H = simtb_pcolor(1:nC, 1:nC, .8*CM);
%     axis square; 
%     axis ij
%     set(gca, 'XTick', [], 'YTick', [], 'CLim', [-1 1])%, 'XColor', [1 1 1], 'YColor', [1 1 1])
%     c = get(gca, 'Children');
%     set(c(find(strcmp(get(c, 'Type'),'line'))), 'Color', 'w');
%     title(sprintf('State %d', ii))
end
cd(Desdir)
save CMAll.mat CMAll

% eT=zeros(nT,nC);
% define the order and time in each state
% Sorder = [1   2  3  4   2]; % state order
% Sdwell = [35 23 40  28  22];

%NOTE: the Sdwell should sum to nT, check here and amend the last partition:
if sum(Sdwell) ~= nT
    Sdwell(end) = nT - sum(Sdwell(1:end-1));
end
Cdwell = cumsum(Sdwell);
Cdwell = [0;Cdwell];
STATE = zeros(1,nT); % state vector
for ii = 1:length(Sorder)
    sIND = Cdwell(ii)+1:Cdwell(ii+1);
    % events related to each module
    rng('default')
    e = rand(length(sIND),nE) < pState;
    e = e.*sign(rand(length(sIND), nE)-0.5);
    for cc = 1:nC
        eTT(sIND,cc) = sign(ModMem(cc,Sorder(ii)))*e(:,abs(ModMem(cc,Sorder(ii))));
    end
    STATE(sIND) = Sorder(ii);
end
cd(Desdir)
save Sdwell.mat Sdwell
save Sorder.mat Sorder
save STATE.mat STATE
P(1) = 6;     % delay of response (relative to onset)
P(2) = 15;    % delay of undershoot (relative to onset)
P(3) = 1;     % dispersion of response
P(4) = 1;     % dispersion of undershoot
P(5) = 3;     % ratio of response to undershoot
P(6) = 0;     % onset (seconds)
P(7) = 32;    % length of kernel (seconds)
TCT  = zeros(nT,nC);
for cc = 1:nC
    TCT(:,cc) = simtb_TCsource(eTT(:,cc), TR, 1, P); % all use same HRF
%     TC(:,cc) = simtb_TCsource(eT(:,cc), TR, 1); % different HRFs
end
%% Create the event time courses
for sub=1:Simsubs
% random aspects (different for each component)
eTR = rand(nT, nC) < pU;
eTR = eTR.*sign(rand(nT, nC)-0.5);
eTR = eTR*aU;

% event time series are stored in eT
%% Convolve event TCs for rest
% [tc, MDESC, P, PDESC] = simtb_TCsource(eT(:,1), TR, 1);

TCR  = zeros(nT,nC);
for cc = 1:nC
    TCR(:,cc) = simtb_TCsource(eTR(:,cc), TR, 1, P); % all use same HRF
%     TC(:,cc) = simtb_TCsource(eT(:,cc), TR, 1); % different HRFs
end

% Add a little gaussian noise
TC = TCR*3+TCT+NoiseSD*randn(nT,nC);
SimsubTCs(:,sub,:)=TC;
% %% Figure to display the states, TCs, and correlation matrices for each partition
% F=figure('color','w','Name', 'sim_TCs_CorrMatrices'); 
% subplot(4, length(Sorder), 1:length(Sorder))
% plot((0:nT-1)*TR, STATE , 'k', 'Linewidth', 1); axis tight; box off
% ylabel('State')
% set(gca, 'YTick', 1:nStates, 'XTick', Cdwell*TR, 'TickDir', 'out', 'Layer', 'Bottom'); grid on
% 
% subplot(4, length(Sorder), length(Sorder)+1:length(Sorder)*2)
% plot((0:nT-1)*TR, TC, 'LineWidth',0.75);
% xlabel('Time (s)')
% ylabel('Amplitude')
% set(gca, 'TickDir', 'out', 'XTick', Cdwell*TR, 'Xgrid', 'on'); 
% axis tight; box off
% 
% for ii = 1:length(Sorder)
%     subplot(4,length(Sorder),length(Sorder)*3+ii)
%     sIND = Cdwell(ii)+1:Cdwell(ii+1);
%     temp = corr(TC(sIND,:));
%     H = simtb_pcolor(1:nC, 1:nC, temp);
%     axis square; axis ij 
%     set(gca, 'XTick', [], 'YTick', [], 'CLim', [-1 1])%, 'XColor', [1 1 1], 'YColor', [1 1 1])
%     c = get(gca, 'Children');
%     set(c(find(strcmp(get(c, 'Type'),'line'))), 'Color', 'w');        
%     text(1.5,-2,sprintf('Partition %d\nState %d', ii, Sorder(ii)), 'Fontsize', 8);
% end
end
cd(Desdir)
save SimsubTCs.mat SimsubTCs