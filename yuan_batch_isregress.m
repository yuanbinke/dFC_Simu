clear all;clc;
rootdir='I:\ISDCC\Newsimtb_0918';
% atlasdir='G:\Narratives\Simony_2016_data\Atlas';
method='L1';
TR=2;
N_time=148;
allpair = 0; parallel = 0;
N_sub=100;
N_roi=[10];
Nwin = N_time;
atlastype={'Simu10.nii'};
datatype={'FixHRF_pu01_au05_NoiseSD0.1','FixHRF_pu01_au05_NoiseSD0.3','FixHRF_pu01_au05_NoiseSD0.6','NoiseSD_FIXHRF_pu01_au05_Poisson2345_NoiseSD0.1','NoiseSD_FIXHRF_pu01_au05_Poisson2345_NoiseSD0.3','NoiseSD_FIXHRF_pu01_au05_Poisson2345_NoiseSD0.6'};

%% dynamic FC
for a=1:length(atlastype)
    for d=4:length(datatype)
        tcdir=fullfile(['I:\ISDCC\Newsimtb_0918\dFC_nozscore_1TR_simu\HMM_HSMM_mhsmm' filesep datatype{d} filesep atlastype{a}(1:end-4) filesep 'tc']);
        cd(tcdir)
        resultdir=tcdir;
        tcfn=fullfile(tcdir,'SimsubTCs.mat');
%         load(tcfn)
%         data=eval('total_tc');
%         clear total_tc
%         data=reshape(data,[N_time+20,N_sub,N_roi(a)]);
        data=importdata(tcfn);
        
        SimsubTCs_isres=[];
       for s=1:N_sub
            
            %is_dcc
            subtc=squeeze(data(:,s,:));%time * ROI
            subR=data;
            subR(:,s,:)=[];
            subR2=squeeze(mean(subR,2));
%             subtc2=[subtc,squeeze(mean(subR,2))];
%             subtc2=subtc2(:,:);
             for ssi=1:size(subtc,2)
            [~,res(:,ssi),~,~, ~, ~, ~] = y_regress_ss(subtc(:,ssi),subR2(:,ssi),1,'T');
             end
            SimsubTCs_isres(:,s,:)=subtc-res;
       end
       save SimsubTCs_isres.mat SimsubTCs_isres
    end
end