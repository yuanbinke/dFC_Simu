clear all;clc;
rootdir='I:\ISDCC\Newsimtb_0918';
% atlasdir='G:\Narratives\Simony_2016_data\Atlas';
method='L1';
TR=2;
wsize=1;
N_time=5704;
N_sub=100;
N_roi=[10];
Nwin = N_time - wsize+1;
atlastype={'Simu10.nii'};
datatype={'NoiseSD_FIXHRF_pu01_au05_Poisson2345_NoiseSD0.1','NoiseSD_FIXHRF_pu01_au05_Poisson2345_NoiseSD0.3','NoiseSD_FIXHRF_pu01_au05_Poisson2345_NoiseSD0.6'};

MTDwsize=5;
% for a=1:length(N_roi)
%     Vtem = spm_vol([atlasdir filesep atlastype{a}]);
%     [Ytem, ~] = spm_read_vols(Vtem);
% %     Ytem(isnan(Ytem)) = 0;
% %     Ytem=round(Ytem);
%     RoiIndex=1:N_roi(a);
%     MNI_coord = cell(length(RoiIndex),1);
%     for j = 1:length(RoiIndex)
%         Region = RoiIndex(j);
%         ind = find(Region == Ytem(:));
% 
%             if ~isempty(ind)
%         [I,J,K] = ind2sub(size(Ytem),ind);
%         XYZ = [I J K]';
%         XYZ(4,:) = 1;
%         MNI_coord{j,1} = XYZ;
%             else
%                 error (['There are no voxels in ROI' blanks(1) num2str(RoiIndex(j)) ', please specify ROIs again']);
%             end
%     end
%     for d=1:length(datatype)
%         tcdir=fullfile([rootdir filesep 'dFC_nozscore_1TR_Intact' filesep datatype{d} filesep atlastype{a}(1:end-4) filesep 'tc']);mkdir(tcdir)
%         datapath=fullfile(rootdir,datatype{d});
%         cd(datapath)
%         sublist=dir('In*');
%         total_z_data=[];
%         for s=1:length(sublist)
% %             [data,header]=y_Read(fullfile(datapath,sublist(s).name,'Filtered_4DVolume.nii'));
% %             sublist(s).name
% %             [M,N,O,T]=size(data);
% %             data=reshape(data,[M*N*O,T]);
%             %         %%voxel z score
%             % %         z_data=zscore(data')';
%             %         z_data=data;
%             %         clear data
%             % %         [Datlas,~]=y_Read(fullfile(roidir,Atlaslist(s-2).name));
%             %         for i = 1:N_roi(a)
%             %             roi=zeros(61,73,61);
%             %             roi(temmask==i)=i;
%             %             roi=roi(:);
%             %             idx=find(roi);
%             %             z_data_roi=z_data(idx,:);
%             %             roi_tc(i,:)=mean(z_data_roi,1);
%             %         end
%             %         %%roi z_score
%             % %         z_roi_tc=zscore(roi_tc');
%             % %         total_z_data=[total_z_data;z_roi_tc];
% 
% 
%             fprintf('Extracting time series for %s\n', sublist(s).name);
%             File_filter='';
%             cd ([datapath filesep sublist(s).name])
%             File_name = spm_select('List',pwd, ['^' File_filter '.*\.img$']);
%             if isempty(File_name)
%                 File_name = spm_select('List',pwd, ['^' File_filter '.*\.nii$']);
%             end
% 
%             Vin = spm_vol(File_name);
%             MTC = zeros(size(Vin,1),length(RoiIndex));
% 
%             for j = 1:length(RoiIndex)
%                 VY = spm_get_data(Vin,MNI_coord{j,1});
%                 MTC(:,j) = mean(VY,2);
%             end
%             MTC(isnan(MTC))=0;
%             total_z_data=[total_z_data;MTC];
%         end
%         tcname=fullfile(tcdir,'total_tc.txt');
%         save(tcname,'total_z_data','-ASCII','-DOUBLE','-TABS')
%     end
% end
% %% static FC
% for a=1:length(N_roi)
%     for d=1:length(datatype)
%         tcdir=fullfile([rootdir filesep 'dFC_nozscore_1TR' filesep datatype{d} filesep atlastype{a}(1:end-4) filesep 'tc']);
%         cd(tcdir)
%         resultdir=fullfile([rootdir filesep 'dFC_nozscore_1TR' filesep datatype{d} filesep atlastype{a}(1:end-4) filesep 'sFC']);mkdir(resultdir)
%         tcfn=fullfile(tcdir,'total_tc.txt');
%         load(tcfn)
%         data=eval('total_tc');
%         clear total_tc
%         %
%         data=reshape(data,[N_time+20,N_sub,N_roi(a)]);
%         dFC_result=[];
%         for s=1:N_sub
%             s
%             subtc=squeeze(data(10:290,s,:));%time * ROI
%             corr_matrix=corr(subtc,subtc);
%             z(:,:,s)=0.5*log((1+corr_matrix)./(1-corr_matrix));
%         end%s
%         ave_z=squeeze(mean(z,3));
%         ave_z(isinf(ave_z))=0;
%         imagesc(ave_z)
%         colorbar
%         colormap jet
%         figurename=fullfile(resultdir,'ave_corr_matrix.jpg');
%         saveas(gcf,figurename);
%         close(gcf)
%         cd(resultdir)
%         save('ave_z.mat','ave_z')
%         clear z
%     end
% end

%% dynamic FC
for a=1:length(atlastype)
    for d=1:length(datatype)
        tcdir=fullfile([rootdir filesep 'dFC_nozscore_1TR_simu' filesep datatype{d} filesep atlastype{a}(1:end-4) filesep 'tc']);
        cd(tcdir)
        resultdir=fullfile([rootdir filesep 'dFC_nozscore_1TR_simu' filesep datatype{d} filesep atlastype{a}(1:end-4) filesep 'ISMTD_Z_wsize' num2str(MTDwsize)]);mkdir(resultdir)
        tcfn=fullfile(tcdir,'SimsubTCs.mat');
%         load(tcfn)
%         data=eval('total_tc');
%         clear total_tc
%         data=reshape(data,[N_time+20,N_sub,N_roi(a)]);
        data=importdata(tcfn);
        
        dFC_result=[];
        for s=1:N_sub
            
            %is_dcc
            subtc=squeeze(data(:,s,:));%time * ROI
            subR=data;
            subR(:,s,:)=[];
            subtc2=[subtc,squeeze(mean(subR,2))];
            subtc2=subtc2(:,:);
             subtc2Z=zscore(subtc2);
            %[tmp_dFC]=pp_ReHo_dALFF_dFC_gift(subtc,method,TR,wsize);%trme * ROI paris, 2D. r*(r-1)/2
%             [~,Ct2,~,~] = DCC_X(subtc2,allpair, parallel);
             Ct2 = coupling(subtc2Z,MTDwsize);
            % extract the upper right ISDCC values
            ISCt2=Ct2(1:N_roi(a),N_roi(a)+1:N_roi(a)*2,:);
%             ISCt2=zeros(N_roi(a),N_roi(a),N_time);
%             for i = 1:N_roi(a)
%                 for j = N_roi(a)+1:N_roi(a)*2
%                     ISCt2(i,j-N_roi,:) = Ct2(i,j,:);
%                 end
%             end
            % moving average DCC with window length
            tmp_dFC_DCCX=zeros(Nwin,N_roi(a)*N_roi(a));
            for iw=1:Nwin
                tmpr=ISCt2(:,:,iw);
                tmp_dFC_DCCX(iw,:)=mat2vec_Asym(tmpr);
            end
            tmp_dFC=tmp_dFC_DCCX(1:Nwin-1,:);
            DEV = std(tmp_dFC, [], 2);%STD OF NODE
            [xmax, imax, xmin, imin] = icatb_extrema(DEV);%local maxima in FC variance
            pIND = sort(imax);%?
            k1_peaks(s) = length(pIND);%?
            SP{s,1} = tmp_dFC(pIND, :);%Subsampling
            dFC_result=[dFC_result;tmp_dFC];
        end%s
        cd(resultdir)
        save('SP.mat','SP','-v7.3')
        
        save('dFC_result.mat','dFC_result','-v7.3')
    end
end

%% clustering
for a=1:length(N_roi)
    for d=1:length(datatype)
        resultdir=fullfile([rootdir filesep 'dFC_nozscore_1TR_simu' filesep datatype{d} filesep atlastype{a}(1:end-4) filesep 'ISMTD_Z_wsize' num2str(MTDwsize)]);
        resultdir2=fullfile([rootdir filesep 'dFC_nozscore_1TR_simu' filesep datatype{d} filesep atlastype{a}(1:end-4) filesep 'kmeans_elbow_ISMTD_Z_wsize' num2str(MTDwsize)]);mkdir(resultdir2)
        kmeans_max_iter = 150;
        dmethod = 'city';
        kmeans_num_replicates = 5;
        num_tests_est_clusters = 10;
        cd(resultdir)
        load SP
        load dFC_result
        %% Cluster
        SPflat = cell2mat(SP);
        clear SP;
        cd(resultdir2)
        cluster_estimate_results = icatb_optimal_clusters(SPflat, min([max(size(SPflat)), 10]), 'method', 'elbow', 'cluster_opts', {'Replicates', kmeans_num_replicates, 'Distance', dmethod, ...
            'MaxIter', kmeans_max_iter}, 'num_tests', num_tests_est_clusters, 'display', 1);
        
        num_clusters = 0;
        for Ntests = 1:length(cluster_estimate_results)
            num_clusters = num_clusters + (cluster_estimate_results{Ntests}.K(1));
        end
        num_clusters = ceil(num_clusters/length(cluster_estimate_results));
        disp(['Number of estimated clusters used in dFNC standard analysis is mean of all tests: ', num2str(num_clusters)]);
        fprintf('\n');
        
        [IDXp, Cp, SUMDp, Dp] = kmeans(SPflat, num_clusters, 'distance', dmethod, 'Replicates', kmeans_num_replicates, 'MaxIter', kmeans_max_iter, 'Display', 'iter', 'empty', 'drop');%gift 4.0b
        
        [IDXall, Call, SUMDall, Dall] = kmeans(dFC_result, num_clusters, 'distance', dmethod, 'Replicates', 1, 'Display', 'iter', 'MaxIter', kmeans_max_iter, ...
            'empty', 'drop', 'Start', Cp);
        
        Tmpmin=zeros(size(Call,1),1);
        Tmpmax=zeros(size(Call,1),1);
        for i=1:size(Call,1)
            tmp_state=sf_vec2mat_Asy(N_roi(a),Call(i,:));
%             tmp_state=tmp_state+tmp_state';
            Tmpmin(i)=min(min(tmp_state));
            Tmpmax(i)=max(max(tmp_state));
        end
        for i=1:size(Call,1)
            tmp_state=sf_vec2mat_Asy(N_roi(a),Call(i,:));
%             tmp_state=tmp_state+tmp_state';
            figure
            imagesc(tmp_state)
            colormap jet
            colorbar
            caxis([min(Tmpmin), max(Tmpmax)]);
            title(['state0' num2str(i)])
            figurename=fullfile(resultdir2,['state0' num2str(i) '.jpg'] );
            saveas(gcf,figurename)
            close(gcf)
            figurename2=strcat('state_0', num2str(i), '.mat') ;
            cd(resultdir2)
            save(figurename2,'tmp_state')
        end
        %
        cd(resultdir2)
        save('IDXall.mat','IDXall')
        save('Call.mat','Call');
        save('SUMDall.mat','SUMDall');
        save('Dall.mat','Dall');
        save('cluster_estimate_results.mat','cluster_estimate_results');
        % %% time parameters
        % kmeansdir=fullfile(rootdir,'kmeans_elbow',datatype{1});
        %         cd(resultdir2)
        %         load('IDXall.mat');
        labels=IDXall;
        %% calulate time
        % TR=2;
        K=max(labels);
        T=length(labels);
        T2=T/N_sub;%number of sliding windows
        for s=1:N_sub
            label_sub=labels(((s-1)*T2+1:s*T2),:);
            dwell_time(s,:)=sf_dwell_time(label_sub,K,TR);
            average_dwell_time(s,:)=sf_ave_dwell_time(label_sub,K,TR);
            transitions_to_state(s,:)=sf_trans_to_state(label_sub,K);
            state_to_state(s,:,:)=sf_state_to_state(label_sub,K);
        end
        cd(resultdir2)
        save dwell_time.mat dwell_time
        save average_dwell_time.mat average_dwell_time
        save transitions_to_state.mat transitions_to_state
        save state_to_state.mat state_to_state
        clear state_to_state transitions_to_state average_dwell_time dwell_time
    end
end
