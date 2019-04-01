
clear all
close all

% loading all data from 3 blocks to plot ERPs for different bands...

% First block; reaching a target point from home: 5cm. each trial: one forward & one backward
% 3 trials for 0degree.
% 2 trials for 45 degree. (toward to body)
% 3 trials for 90 degree. (toward to body)
% It seems chs 63, 95, 109 were bad chs.
%
% Second block; reaching a target point from home: 7.5cm. each trial: one forward & one backward
% 3 trials for 0degree.
% 3 trials for 45 degree. (toward to body)
% 3 trials for 90 degree. (toward to body)
% It seems chs 63, 95, 109 were bad chs.
%
% Third block; reaching a target point from home: 7.5cm. each trial: one forward & one backward
% 3 trials for 0degree.
% 3 trials for 45 degree. (toward to body)
% 3 trials for 90 degree. (toward to body)
% It seems chs 63, 95, 109 were bad chs.

%% loading and breaking raw ECoG data into trials
% first block
load('..\20190321-113426\20190321-113426-001.mat');
timestamp_B1=anin;
ECoG_B1=lfp;
clear anin lfp
% second block
load('..\20190321-113426\20190321-113426-002.mat');
timestamp_B2=anin;
ECoG_B2=lfp;
clear anin lfp
% third block
load('..\20190321-113426\20190321-113426-003.mat');
timestamp_B3=anin;
ECoG_B3=lfp;
clear anin lfp

if 0
    % visualizing check-up; blackrock data; not compatible with grid chs number
    data=ECoG_B3;
    t1=5001;
    t2=6000;
    figure(1);
    for i=1:64
        subplot(8,8,i)
        plot(data(t1:t2,i))
        title(['Ch',num2str(i)])
        
    end
    figure(2);
    for i=65:128
        subplot(8,8,i-64)
        plot(data(t1:t2,i))
        title(['Ch',num2str(i)])
        
    end
    
    plot(data(1001:10000,63))
    
end

% Excluding bad channels
BadChs=[63, 95];
All_Index=ones(size(ECoG_B1,2),1);
All_Index(BadChs')=0;
Selected_Chs=logical(All_Index(1:size(ECoG_B1,2),1));
ECoG_B1(:,BadChs)=0;
ECoG_B2(:,BadChs)=0;
ECoG_B3(:,BadChs)=0;

% z score
ECoG_B1_1=zscore(ECoG_B1)';
ECoG_B2_1=zscore(ECoG_B2)';
ECoG_B3_1=zscore(ECoG_B3)';

% option median or mean:
% F1_Mean =mean(F1_ECoGdata(Selected_Chs,:),1);
% F1_ECoGdata= F1_ECoGdata-repmat(F1_Mean,256,1);

ECoG_B1_Med=median(ECoG_B1_1,1);
ECoG_B1_2= ECoG_B1_1-repmat(ECoG_B1_Med,size(ECoG_B1,2),1);
ECoG_B1_2=ECoG_B1_2';

ECoG_B2_Med=median(ECoG_B2_1,1);
ECoG_B2_2= ECoG_B2_1-repmat(ECoG_B2_Med,size(ECoG_B2,2),1);
ECoG_B2_2=ECoG_B2_2';

ECoG_B3_Med=median(ECoG_B3_1,1);
ECoG_B3_2= ECoG_B3_1-repmat(ECoG_B3_Med,size(ECoG_B3,2),1);
ECoG_B3_2=ECoG_B3_2';

ECoG_B1_2(:,BadChs)=0;
ECoG_B2_2(:,BadChs)=0;
ECoG_B3_2(:,BadChs)=0;

% pulling out Filtered data
FilterBands={
    [.5, 4]% delta
    [4,8]% theta
    [8,13]% alpha
    [13,30]% beta
    [30,58] % gamma 1
    [70,110]}; % gamma 2

% Cal ERPs for each band for each angle; for all blocks
% selecting window for ERP in each band; sample points
ERPwindowPerBand={
    [3000, 3000]% delta
    [1000,1000]% theta
    [500,500]% alpha
    [500,500]% beta
    [500,500] % gamma 1
    [500,500]}; % gamma 2;

%index of timestamps for all blocks
threshhold=1000;
Stim_B1=find(diff(timestamp_B1>threshhold)==1);
index_StartReaching_B1=Stim_B1(1:2:end);
index_BackingHome_B1=Stim_B1(2:2:end);

Stim_B2=find(diff(timestamp_B2>threshhold)==1);
index_StartReaching_B2=Stim_B2(1:2:end);
index_BackingHome_B2=Stim_B2(2:2:end);

Stim_B3=find(diff(timestamp_B3>threshhold)==1);
index_StartReaching_B3=Stim_B3(1:2:end);
index_BackingHome_B3=Stim_B3(2:2:end);

ch_layout = [
    91	84	67	90	70	79	88	69	92	83	65	89	87	86	94	82
    66	93	78	95	76	75	85	73	68	80	74	72	96	71	77	81
    60	37	42	50	56	54	49	40	43	35	45	63	47	46	58	55
    53	57	33	48	39	51	41	34	64	52	62	38	36	44	61	59
    8	26	29	28	9	5	13	20	11	23	16	22	27	4	3	31
    7	21	15	24	25	1	2	32	14	12	30	19	18	17	6	10
    110	125	111	115	103	117	100	123	113	119	118	98	101	105	116	99
    107	112	97	128	121	124	108	109	127	126	106	122	114	120	104	102];

ch_layout=ch_layout';
ch_layout=ch_layout(:);

bandname={'delta','theta','alpha','beta','gamma1','gamma2'};

for band=1 %:6 % including raw(?!)
    
    [b,a]=butter(3,FilterBands{band}/(Fs/2));
    ECoG_B1_Filtered=filtfilt(b,a,ECoG_B1_2);
    ECoG_B2_Filtered=filtfilt(b,a,ECoG_B2_2);
    ECoG_B3_Filtered=filtfilt(b,a,ECoG_B3_2);
    
%     ECoG_B1_Filtered=abs(hilbert(ECoG_B1_Filtered));
%     ECoG_B2_Filtered=abs(hilbert(ECoG_B2_Filtered));
%     ECoG_B3_Filtered=abs(hilbert(ECoG_B3_Filtered));
    
    ERP_SRAll_0=zeros(ERPwindowPerBand{band}(1)+ERPwindowPerBand{band}(2)+1,size(ECoG_B1,2));
    ERP_BHAll_0=zeros(ERPwindowPerBand{band}(1)+ERPwindowPerBand{band}(2)+1,size(ECoG_B1,2));
    ERP_SRAll_45=zeros(ERPwindowPerBand{band}(1)+ERPwindowPerBand{band}(2)+1,size(ECoG_B1,2));
    ERP_BHAll_45=zeros(ERPwindowPerBand{band}(1)+ERPwindowPerBand{band}(2)+1,size(ECoG_B1,2));
    ERP_SRAll_90=zeros(ERPwindowPerBand{band}(1)+ERPwindowPerBand{band}(2)+1,size(ECoG_B1,2));
    ERP_BHAll_90=zeros(ERPwindowPerBand{band}(1)+ERPwindowPerBand{band}(2)+1,size(ECoG_B1,2));
    
    for trial=1:3
        idx_SR_B1=index_StartReaching_B1(trial)-ERPwindowPerBand{band}(1):index_StartReaching_B1(trial)+ERPwindowPerBand{band}(2);
        idx_SR_B2=index_StartReaching_B2(trial)-ERPwindowPerBand{band}(1):index_StartReaching_B2(trial)+ERPwindowPerBand{band}(2);
        idx_SR_B3=index_StartReaching_B3(trial)-ERPwindowPerBand{band}(1):index_StartReaching_B3(trial)+ERPwindowPerBand{band}(2);
        ERP_SR_0=ECoG_B1_Filtered(idx_SR_B1,:)+ECoG_B2_Filtered(idx_SR_B2,:)+ECoG_B3_Filtered(idx_SR_B3,:);
        ERP_SRAll_0=ERP_SRAll_0+ERP_SR_0;
        
        idx_BH_B1=index_BackingHome_B1(trial)-ERPwindowPerBand{band}(1):index_BackingHome_B1(trial)+ERPwindowPerBand{band}(2);
        idx_BH_B2=index_BackingHome_B2(trial)-ERPwindowPerBand{band}(1):index_BackingHome_B2(trial)+ERPwindowPerBand{band}(2);
        idx_BH_B3=index_BackingHome_B3(trial)-ERPwindowPerBand{band}(1):index_BackingHome_B3(trial)+ERPwindowPerBand{band}(2);
        ERP_BH_0=ECoG_B1_Filtered(idx_BH_B1,:)+ECoG_B2_Filtered(idx_BH_B2,:)+ECoG_B3_Filtered(idx_BH_B3,:);
        ERP_BHAll_0=ERP_BHAll_0+ERP_BH_0;
        
        idx_SR_B1=index_StartReaching_B1(trial+5)-ERPwindowPerBand{band}(1):index_StartReaching_B1(trial+5)+ERPwindowPerBand{band}(2);
        idx_SR_B2=index_StartReaching_B2(trial+6)-ERPwindowPerBand{band}(1):index_StartReaching_B2(trial+6)+ERPwindowPerBand{band}(2);
        idx_SR_B3=index_StartReaching_B3(trial+6)-ERPwindowPerBand{band}(1):index_StartReaching_B3(trial+6)+ERPwindowPerBand{band}(2);
        ERP_SR_90=ECoG_B1_Filtered(idx_SR_B1,:)+ECoG_B2_Filtered(idx_SR_B2,:)+ECoG_B3_Filtered(idx_SR_B3,:);
        ERP_SRAll_90=ERP_SRAll_90+ERP_SR_90;
        
        idx_BH_B1=index_BackingHome_B1(trial+5)-ERPwindowPerBand{band}(1):index_BackingHome_B1(trial+5)+ERPwindowPerBand{band}(2);
        idx_BH_B2=index_BackingHome_B2(trial+6)-ERPwindowPerBand{band}(1):index_BackingHome_B2(trial+6)+ERPwindowPerBand{band}(2);
        idx_BH_B3=index_BackingHome_B3(trial+6)-ERPwindowPerBand{band}(1):index_BackingHome_B3(trial+6)+ERPwindowPerBand{band}(2);
        ERP_BH_90=ECoG_B1_Filtered(idx_BH_B1,:)+ECoG_B2_Filtered(idx_BH_B2,:)+ECoG_B3_Filtered(idx_BH_B3,:);
        ERP_BHAll_90=ERP_BHAll_90+ERP_BH_90;
        
        idx_SR_B2=index_StartReaching_B2(trial+3)-ERPwindowPerBand{band}(1):index_StartReaching_B2(trial+3)+ERPwindowPerBand{band}(2);
        idx_SR_B3=index_StartReaching_B3(trial+3)-ERPwindowPerBand{band}(1):index_StartReaching_B3(trial+3)+ERPwindowPerBand{band}(2);
        ERP_SR_45=ECoG_B2_Filtered(idx_SR_B2,:)+ECoG_B3_Filtered(idx_SR_B3,:);
        ERP_SRAll_45=ERP_SRAll_45+ERP_SR_45;
        
        idx_BH_B2=index_BackingHome_B2(trial+3)-ERPwindowPerBand{band}(1):index_BackingHome_B2(trial+3)+ERPwindowPerBand{band}(2);
        idx_BH_B3=index_BackingHome_B3(trial+3)-ERPwindowPerBand{band}(1):index_BackingHome_B3(trial+3)+ERPwindowPerBand{band}(2);
        ERP_BH_45=ECoG_B2_Filtered(idx_BH_B2,:)+ECoG_B3_Filtered(idx_BH_B3,:);
        ERP_BHAll_45=ERP_BHAll_45+ERP_BH_45;
        
    end
    
    for trial=1:2
        
        idx_SR_B1=index_StartReaching_B1(trial+3)-ERPwindowPerBand{band}(1):index_StartReaching_B1(trial+3)+ERPwindowPerBand{band}(2);
        ERP_SRAll_45=ERP_SRAll_45+ECoG_B1_Filtered(idx_SR_B1,:);
        
        idx_BH_B1=index_BackingHome_B1(trial+3)-ERPwindowPerBand{band}(1):index_BackingHome_B1(trial+3)+ERPwindowPerBand{band}(2);
        ERP_BHAll_45=ERP_BHAll_45+ECoG_B1_Filtered(idx_BH_B1,:);
        
    end
    
    ERP_SRAll_0= ERP_SRAll_0/9;
    ERP_BHAll_0= ERP_BHAll_0/9;
    ERP_SRAll_90= ERP_SRAll_90/9;
    ERP_BHAll_90= ERP_BHAll_90/9;
    ERP_SRAll_45= ERP_SRAll_45/8;
    ERP_BHAll_45= ERP_BHAll_45/8;
    
    % for 0 degree
    figure;
    set(gcf, 'Position', [100, 100, 2400, 1200]);
    suptitle(['ERPs for all blocks & 0 degree & Start-Reaching; ',bandname{band},' band'])
    for i=1:128
        subplot(8,16,i)
        plot(ERP_SRAll_0(:,ch_layout(i)))
        xlim([0 ERPwindowPerBand{band}(1)+ERPwindowPerBand{band}(2)]);
        vline(ERPwindowPerBand{band}(1));
        xticks([1 ERPwindowPerBand{band}(1) ERPwindowPerBand{band}(1)+ERPwindowPerBand{band}(2)]);
        xticklabels({num2str(-ERPwindowPerBand{band}(1)),num2str(0),num2str(ERPwindowPerBand{band}(2))});
        title(['Ch',num2str(ch_layout(i))]);
    end
    % Producing high quality iamge and save it
    %HighQualityFigs(['ERPsForAllBlocks_0Degree_SR_',bandname{band}])
    
    figure;
    set(gcf, 'Position', [100, 100, 2400, 1200]);
    suptitle(['ERPs for all blocks & 0 degree & Backing-Home; ',bandname{band},' band'])
    for i=1:128
        subplot(8,16,i)
        plot(ERP_BHAll_0(:,ch_layout(i)))
        xlim([0 ERPwindowPerBand{band}(1)+ERPwindowPerBand{band}(2)]);
        vline(ERPwindowPerBand{band}(1));
        xticks([1 ERPwindowPerBand{band}(1) ERPwindowPerBand{band}(1)+ERPwindowPerBand{band}(2)]);
        xticklabels({num2str(-ERPwindowPerBand{band}(1)),num2str(0),num2str(ERPwindowPerBand{band}(2))});
        title(['Ch',num2str(ch_layout(i))]);
    end
    % Producing high quality iamge and save it
    %HighQualityFigs(['ERPsForAllBlocks_0Degree_BH_',bandname{band}])
    
    figure;
    set(gcf, 'Position', [100, 100, 2400, 1200]);
    suptitle(['ERPs for all blocks & 90 degree & Start-Reaching+Backing-Home; ',bandname{band},' band'])
    for i=1:128
        subplot(8,16,i)
        plot((ERP_SRAll_0(:,ch_layout(i))+ERP_BHAll_0(:,ch_layout(i)))/2)
        xlim([0 ERPwindowPerBand{band}(1)+ERPwindowPerBand{band}(2)]);
        vline(ERPwindowPerBand{band}(1));
        xticks([1 ERPwindowPerBand{band}(1) ERPwindowPerBand{band}(1)+ERPwindowPerBand{band}(2)]);
        xticklabels({num2str(-ERPwindowPerBand{band}(1)),num2str(0),num2str(ERPwindowPerBand{band}(2))});
        title(['Ch',num2str(ch_layout(i))]);
    end
    % Producing high quality iamge and save it
    %HighQualityFigs(['ERPsForAllBlocks_0Degree_SR&BH_',bandname{band}])
    
    % for 45 degree
    figure;
    set(gcf, 'Position', [100, 100, 2400, 1200]);
    suptitle(['ERPs for all blocks & 45 degree & Start-Reaching; ',bandname{band},' band'])
    for i=1:128
        subplot(8,16,i)
        plot(ERP_SRAll_45(:,ch_layout(i)))
        xlim([0 ERPwindowPerBand{band}(1)+ERPwindowPerBand{band}(2)]);
        vline(ERPwindowPerBand{band}(1));
        xticks([1 ERPwindowPerBand{band}(1) ERPwindowPerBand{band}(1)+ERPwindowPerBand{band}(2)]);
        xticklabels({num2str(-ERPwindowPerBand{band}(1)),num2str(0),num2str(ERPwindowPerBand{band}(2))});
        title(['Ch',num2str(ch_layout(i))]);
    end
    % Producing high quality iamge and save it
    %HighQualityFigs(['ERPsForAllBlocks_45Degree_SR_',bandname{band}])
    
    figure;
    set(gcf, 'Position', [100, 100, 2400, 1200]);
    suptitle(['ERPs for all blocks & 45 degree & Backing-Home; ',bandname{band},' band'])
    for i=1:128
        subplot(8,16,i)
        plot(ERP_BHAll_45(:,ch_layout(i)))
        xlim([0 ERPwindowPerBand{band}(1)+ERPwindowPerBand{band}(2)]);
        vline(ERPwindowPerBand{band}(1));
        xticks([1 ERPwindowPerBand{band}(1) ERPwindowPerBand{band}(1)+ERPwindowPerBand{band}(2)]);
        xticklabels({num2str(-ERPwindowPerBand{band}(1)),num2str(0),num2str(ERPwindowPerBand{band}(2))});
        title(['Ch',num2str(ch_layout(i))]);
    end
    % Producing high quality iamge and save it
    %HighQualityFigs(['ERPsForAllBlocks_45Degree_BH_',bandname{band}])
    
    figure;
    set(gcf, 'Position', [100, 100, 2400, 1200]);
    suptitle(['ERPs for all blocks & 45 degree & Start-Reaching+Backing-Home; ',bandname{band},' band'])
    for i=1:128
        subplot(8,16,i)
        plot((ERP_SRAll_45(:,ch_layout(i))+ERP_BHAll_45(:,ch_layout(i)))/2)
        xlim([0 ERPwindowPerBand{band}(1)+ERPwindowPerBand{band}(2)]);
        vline(ERPwindowPerBand{band}(1));
        xticks([1 ERPwindowPerBand{band}(1) ERPwindowPerBand{band}(1)+ERPwindowPerBand{band}(2)]);
        xticklabels({num2str(-ERPwindowPerBand{band}(1)),num2str(0),num2str(ERPwindowPerBand{band}(2))});
        title(['Ch',num2str(ch_layout(i))]);
    end
    % Producing high quality iamge and save it
    %HighQualityFigs(['ERPsForAllBlocks_45Degree_SR&BH_',bandname{band}])
    
    % for 90 degree
    figure;
    set(gcf, 'Position', [100, 100, 2400, 1200]);
    suptitle(['ERPs for all blocks & 0 degree & Start-Reaching; ',bandname{band},' band'])
    for i=1:128
        subplot(8,16,i)
        plot(ERP_SRAll_90(:,ch_layout(i)))
        xlim([0 ERPwindowPerBand{band}(1)+ERPwindowPerBand{band}(2)]);
        vline(ERPwindowPerBand{band}(1));
        xticks([1 ERPwindowPerBand{band}(1) ERPwindowPerBand{band}(1)+ERPwindowPerBand{band}(2)]);
        xticklabels({num2str(-ERPwindowPerBand{band}(1)),num2str(0),num2str(ERPwindowPerBand{band}(2))});
        title(['Ch',num2str(ch_layout(i))]);
    end
    % Producing high quality iamge and save it
    %HighQualityFigs(['ERPsForAllBlocks_90Degree_SR_',bandname{band}])
    
    figure;
    set(gcf, 'Position', [100, 100, 2400, 1200]);
    suptitle(['ERPs for all blocks & 0 degree & Backing-Home; ',bandname{band},' band'])
    for i=1:128
        subplot(8,16,i)
        plot(ERP_BHAll_90(:,ch_layout(i)))
        xlim([0 ERPwindowPerBand{band}(1)+ERPwindowPerBand{band}(2)]);
        vline(ERPwindowPerBand{band}(1));
        xticks([1 ERPwindowPerBand{band}(1) ERPwindowPerBand{band}(1)+ERPwindowPerBand{band}(2)]);
        xticklabels({num2str(-ERPwindowPerBand{band}(1)),num2str(0),num2str(ERPwindowPerBand{band}(2))});
        title(['Ch',num2str(ch_layout(i))]);
    end
    % Producing high quality iamge and save it
    %HighQualityFigs(['ERPsForAllBlocks_90Degree_BH_',bandname{band}])
    
    figure;
    set(gcf, 'Position', [100, 100, 2400, 1200]);
    suptitle(['ERPs for all blocks & 0 degree & Start-Reaching+Backing-Home; ',bandname{band},' band'])
    for i=1:128
        subplot(8,16,i)
        plot((ERP_SRAll_90(:,ch_layout(i))+ERP_BHAll_0(:,ch_layout(i)))/2)
        xlim([0 ERPwindowPerBand{band}(1)+ERPwindowPerBand{band}(2)]);
        vline(ERPwindowPerBand{band}(1));
        xticks([1 ERPwindowPerBand{band}(1) ERPwindowPerBand{band}(1)+ERPwindowPerBand{band}(2)]);
        xticklabels({num2str(-ERPwindowPerBand{band}(1)),num2str(0),num2str(ERPwindowPerBand{band}(2))});
        title(['Ch',num2str(ch_layout(i))]);
    end
    % Producing high quality iamge and save it
    %HighQualityFigs(['ERPsForAllBlocks_90Degree_SR&BH_',bandname{band}])
    
    
    % for all angles 
    figure;
    set(gcf, 'Position', [100, 100, 2400, 1200]);
    suptitle(['ERPs for all blocks for all angles & Start-Reaching; ',bandname{band},' band'])
    for i=1:128
        subplot(8,16,i)
        plot((ERP_SRAll_0(:,ch_layout(i))+ERP_SRAll_45(:,ch_layout(i))+ERP_SRAll_90(:,ch_layout(i)))/3)
        xlim([0 ERPwindowPerBand{band}(1)+ERPwindowPerBand{band}(2)]);
        vline(ERPwindowPerBand{band}(1));
        xticks([1 ERPwindowPerBand{band}(1) ERPwindowPerBand{band}(1)+ERPwindowPerBand{band}(2)]);
        xticklabels({num2str(-ERPwindowPerBand{band}(1)),num2str(0),num2str(ERPwindowPerBand{band}(2))});
        title(['Ch',num2str(ch_layout(i))]);
        
    end
    % Producing high quality iamge and save it
    %HighQualityFigs(['ERPsForAllBlocks&Angles_SR_',bandname{band}])
    
    figure;
    set(gcf, 'Position', [100, 100, 2400, 1200]);
    suptitle(['ERPs for all blocks for all angles & Backing-Home; ',bandname{band},' band'])
    for i=1:128
        subplot(8,16,i)
        plot((ERP_BHAll_0(:,ch_layout(i))+ERP_BHAll_45(:,ch_layout(i))+ERP_BHAll_90(:,ch_layout(i)))/3)
        xlim([0 ERPwindowPerBand{band}(1)+ERPwindowPerBand{band}(2)]);
        vline(ERPwindowPerBand{band}(1));
        xticks([1 ERPwindowPerBand{band}(1) ERPwindowPerBand{band}(1)+ERPwindowPerBand{band}(2)]);
        xticklabels({num2str(-ERPwindowPerBand{band}(1)),num2str(0),num2str(ERPwindowPerBand{band}(2))});
        title(['Ch',num2str(ch_layout(i))]);
        
    end
    % Producing high quality iamge and save it
    %HighQualityFigs(['ERPsForAllBlocks&Angles_BH_',bandname{band}])
    
    figure;
    set(gcf, 'Position', [100, 100, 2400, 1200]);
    suptitle(['ERPs for all blocks for all angles & Start-Reaching+Backing-Home; ',bandname{band},' band'])
    for i=1:128
        subplot(8,16,i)
        plot((ERP_SRAll_0(:,ch_layout(i))+ERP_BHAll_0(:,ch_layout(i))+...
            ERP_SRAll_45(:,ch_layout(i))+ERP_BHAll_45(:,ch_layout(i))+...
            ERP_SRAll_90(:,ch_layout(i))+ERP_BHAll_90(:,ch_layout(i)))/6)
        xlim([0 ERPwindowPerBand{band}(1)+ERPwindowPerBand{band}(2)]);
        vline(ERPwindowPerBand{band}(1));
        xticks([1 ERPwindowPerBand{band}(1) ERPwindowPerBand{band}(1)+ERPwindowPerBand{band}(2)]);
        xticklabels({num2str(-ERPwindowPerBand{band}(1)),num2str(0),num2str(ERPwindowPerBand{band}(2))});
        title(['Ch',num2str(ch_layout(i))]);
        
    end
    % Producing high quality iamge and save it
    %HighQualityFigs(['ERPsForAllBlocks&Angles_SR&BH_',bandname{band}])
    
    
    
    
end








