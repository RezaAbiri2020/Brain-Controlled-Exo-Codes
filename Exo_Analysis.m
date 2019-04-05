
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

% loading and breaking raw ECoG data into trials
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
    [500, 1500]% delta
    [1000,1000]% theta
    [500,500]% alpha
    [500,500]% beta
    [50,100] % gamma 1
    [50,100]}; % gamma 2;

%index of timestamps for all blocks
% E
threshhold=1000;
Stim_B1=find(diff(timestamp_B1>threshhold)==1);
index_E_B1=Stim_B1(1:2:end);
index_F_B1=Stim_B1(2:2:end);

Stim_B2=find(diff(timestamp_B2>threshhold)==1);
index_E_B2=Stim_B2(1:2:end);
index_F_B2=Stim_B2(2:2:end);

Stim_B3=find(diff(timestamp_B3>threshhold)==1);
index_E_B3=Stim_B3(1:2:end);
index_F_B3=Stim_B3(2:2:end);

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

ECoG_B1_3=zeros(size(ECoG_B1_2));
ECoG_B2_3=zeros(size(ECoG_B2_2));
ECoG_B3_3=zeros(size(ECoG_B3_2));

for i=1:128
    ECoG_B1_3(:,i)=ECoG_B1_2(:,ch_layout(i));
    ECoG_B2_3(:,i)=ECoG_B2_2(:,ch_layout(i));
    ECoG_B3_3(:,i)=ECoG_B3_2(:,ch_layout(i));
end 


bandname={'delta','theta','alpha','beta','gamma1','gamma2'};

% for colors adjustment in clims
aclim=-0.4;
bclim=0.4;



for band=2 %:6 % including raw(?!)
    
    [b,a]=butter(3,FilterBands{band}/(Fs/2));
    ECoG_B1_Filtered=filtfilt(b,a,ECoG_B1_3);
    ECoG_B2_Filtered=filtfilt(b,a,ECoG_B2_3);
    ECoG_B3_Filtered=filtfilt(b,a,ECoG_B3_3);
    
    %     ECoG_B1_Filtered=abs(hilbert(ECoG_B1_Filtered));
    %     ECoG_B2_Filtered=abs(hilbert(ECoG_B2_Filtered));
    %     ECoG_B3_Filtered=abs(hilbert(ECoG_B3_Filtered));
    % Generating ERPs for trials in Blocks
    
    for i=1:3 % for number of extension and flexion
        % block 1
        B1(i).E0=ECoG_B1_Filtered(index_E_B1(i)-ERPwindowPerBand{band}(1):index_E_B1(i)+ERPwindowPerBand{band}(2),:);
        B1(i).E90=ECoG_B1_Filtered(index_E_B1(i+5)-ERPwindowPerBand{band}(1):index_E_B1(i+5)+ERPwindowPerBand{band}(2),:);
        
        B1(i).F0=ECoG_B1_Filtered(index_F_B1(i)-ERPwindowPerBand{band}(1):index_F_B1(i)+ERPwindowPerBand{band}(2),:);
        B1(i).F90=ECoG_B1_Filtered(index_F_B1(i+5)-ERPwindowPerBand{band}(1):index_F_B1(i+5)+ERPwindowPerBand{band}(2),:);
        
        %Block 2
        B2(i).E0=ECoG_B2_Filtered(index_E_B2(i)-ERPwindowPerBand{band}(1):index_E_B2(i)+ERPwindowPerBand{band}(2),:);
        B2(i).E45=ECoG_B2_Filtered(index_E_B2(i+3)-ERPwindowPerBand{band}(1):index_E_B2(i+3)+ERPwindowPerBand{band}(2),:);
        B2(i).E90=ECoG_B2_Filtered(index_E_B2(i+6)-ERPwindowPerBand{band}(1):index_E_B2(i+6)+ERPwindowPerBand{band}(2),:);
        
        B2(i).F0=ECoG_B2_Filtered(index_F_B2(i)-ERPwindowPerBand{band}(1):index_F_B2(i)+ERPwindowPerBand{band}(2),:);
        B2(i).F45=ECoG_B2_Filtered(index_F_B2(i+3)-ERPwindowPerBand{band}(1):index_F_B2(i+3)+ERPwindowPerBand{band}(2),:);
        B2(i).F90=ECoG_B2_Filtered(index_F_B2(i+6)-ERPwindowPerBand{band}(1):index_F_B2(i+6)+ERPwindowPerBand{band}(2),:);
        
        %Block 3
        B3(i).E0=ECoG_B3_Filtered(index_E_B3(i)-ERPwindowPerBand{band}(1):index_E_B3(i)+ERPwindowPerBand{band}(2),:);
        B3(i).E45=ECoG_B3_Filtered(index_E_B3(i+3)-ERPwindowPerBand{band}(1):index_E_B3(i+3)+ERPwindowPerBand{band}(2),:);
        B3(i).E90=ECoG_B3_Filtered(index_E_B3(i+6)-ERPwindowPerBand{band}(1):index_E_B3(i+6)+ERPwindowPerBand{band}(2),:);
        
        B3(i).F0=ECoG_B3_Filtered(index_F_B3(i)-ERPwindowPerBand{band}(1):index_F_B3(i)+ERPwindowPerBand{band}(2),:);
        B3(i).F45=ECoG_B3_Filtered(index_F_B3(i+3)-ERPwindowPerBand{band}(1):index_F_B3(i+3)+ERPwindowPerBand{band}(2),:);
        B3(i).F90=ECoG_B3_Filtered(index_F_B3(i+6)-ERPwindowPerBand{band}(1):index_F_B3(i+6)+ERPwindowPerBand{band}(2),:);
             
    end
    
    for i=1:2
        %block 2
        B1(i).E45=ECoG_B1_Filtered(index_E_B1(i+3)-ERPwindowPerBand{band}(1):index_E_B1(i+3)+ERPwindowPerBand{band}(2),:);
        B1(i).F45=ECoG_B1_Filtered(index_F_B1(i+3)-ERPwindowPerBand{band}(1):index_F_B1(i+3)+ERPwindowPerBand{band}(2),:);
        
    end
    
    % inter-block corr for each channel for only extension mode at 0 degree
    
    % block 1 /2 /3
    for i=1:128
        interB_B1_E0(i)=(corr(B1(1).E0(:,i),B1(2).E0(:,i))+corr(B1(1).E0(:,i),B1(3).E0(:,i))+corr(B1(2).E0(:,i),B1(3).E0(:,i)))/3;
        interB_B2_E0(i)=(corr(B2(1).E0(:,i),B2(2).E0(:,i))+corr(B2(1).E0(:,i),B2(3).E0(:,i))+corr(B2(2).E0(:,i),B2(3).E0(:,i)))/3;
        interB_B3_E0(i)=(corr(B3(1).E0(:,i),B3(2).E0(:,i))+corr(B3(1).E0(:,i),B3(3).E0(:,i))+corr(B3(2).E0(:,i),B3(3).E0(:,i)))/3;
    end
    
    % intra-block corr for each channel for only extension mode at 0 degree
    for i=1:128
        
        intraBAll=0;
        for j=1:3
            for jj=1:3
                intraB=(corr(B1(j).E0(:,i),B2(jj).E0(:,i))+corr(B1(j).E0(:,i),B3(jj).E0(:,i))+corr(B2(j).E0(:,i),B3(jj).E0(:,i)))/3;
                intraBAll=intraBAll+intraB;
                
            end
            
        end
        intraB_E0(i)=intraBAll/9;
    end
    
    figure;
    set(gcf, 'Position', [100, 100, 2400, 600]);
    suptitle(['Correlation for each channel; ',bandname{band},' band'])
    subplot(4,1,1)
    clims = [aclim bclim];
    imagesc(interB_B1_E0,clims)
    colorbar
    title('interB-B1-E0');
    
    subplot(4,1,2)
    clims = [aclim bclim];
    imagesc(interB_B2_E0,clims)
    colorbar
    title('interB-B2-E0');
    
    subplot(4,1,3)
    clims = [aclim bclim];
    imagesc(interB_B3_E0,clims)
    colorbar
    title('interB-B3-E0');
    
    subplot(4,1,4)
    clims = [aclim bclim];
    imagesc(intraB_E0,clims)
    colorbar
    title('intraB-E0');
    
   
    figure;
    set(gcf, 'Position', [100, 100, 2400, 1200]);
    suptitle(['Correlation for each channel; ',bandname{band},' band'])
    interB_B1_E0_1=reshape(interB_B1_E0,16,8);
    interB_B1_E0_1=interB_B1_E0_1';
    clims = [aclim bclim];
    subplot(2,2,1)
    imagesc(interB_B1_E0_1,clims)
    colorbar
    title('interB-B1-E0');
    interB_B2_E0_1=reshape(interB_B2_E0,16,8);
    interB_B2_E0_1=interB_B2_E0_1';
    clims = [aclim bclim];
    subplot(2,2,2)
    imagesc(interB_B2_E0_1,clims)
    colorbar
    title('interB-B2-E0');
    interB_B3_E0_1=reshape(interB_B3_E0,16,8);
    interB_B3_E0_1=interB_B3_E0_1';
    clims = [aclim bclim];
    subplot(2,2,3)
    imagesc(interB_B3_E0_1,clims)
    colorbar
    title('interB-B3-E0');
    intraB_E0_1=reshape(intraB_E0,16,8);
    intraB_E0_1=intraB_E0_1';
    clims = [aclim bclim];
    subplot(2,2,4)
    imagesc(intraB_E0_1,clims)
    colorbar
    title('intraB-E0');
    
    
    
    % inter-block corr for each channel for only flexion mode at 0 degree
    % block 1 /2 /3
    for i=1:128
        interB_B1_F0(i)=(corr(B1(1).F0(:,i),B1(2).F0(:,i))+corr(B1(1).F0(:,i),B1(3).F0(:,i))+corr(B1(2).F0(:,i),B1(3).F0(:,i)))/3;
        interB_B2_F0(i)=(corr(B2(1).F0(:,i),B2(2).F0(:,i))+corr(B2(1).F0(:,i),B2(3).F0(:,i))+corr(B2(2).F0(:,i),B2(3).F0(:,i)))/3;
        interB_B3_F0(i)=(corr(B3(1).F0(:,i),B3(2).F0(:,i))+corr(B3(1).F0(:,i),B3(3).F0(:,i))+corr(B3(2).F0(:,i),B3(3).F0(:,i)))/3;
    end
    
    
    % intra-block corr for each channel for only flexion mode at 0 degree
    for i=1:128
        
        intraBAll=0;
        for j=1:3
            for jj=1:3
                intraB=(corr(B1(j).F0(:,i),B2(jj).F0(:,i))+corr(B1(j).F0(:,i),B3(jj).F0(:,i))+corr(B2(j).F0(:,i),B3(jj).F0(:,i)))/3;
                intraBAll=intraBAll+intraB;
                
            end
            
        end
        intraB_F0(i)=intraBAll/9;
    end
   
    figure;
    set(gcf, 'Position', [100, 100, 2400, 600]);
    suptitle(['Correlation for each channel; ',bandname{band},' band'])
    subplot(4,1,1)
    clims = [aclim bclim];
    imagesc(interB_B1_F0,clims)
    colorbar
    title('interB-B1-F0');
    
    subplot(4,1,2)
    clims = [aclim bclim];
    imagesc(interB_B2_F0,clims)
    colorbar
    title('interB-B2-F0');
    
    subplot(4,1,3)
    clims = [aclim bclim];
    imagesc(interB_B3_F0,clims)
    colorbar
    title('interB-B3-F0');
    
    subplot(4,1,4)
    clims = [aclim bclim];
    imagesc(intraB_F0,clims)
    colorbar
    title('intraB-F0');
    
    figure;
    set(gcf, 'Position', [100, 100, 2400, 1200]);
    suptitle(['Correlation for each channel; ',bandname{band},' band'])
    interB_B1_F0_1=reshape(interB_B1_F0,16,8);
    interB_B1_F0_1=interB_B1_F0_1';
    clims = [aclim bclim];
    subplot(2,2,1)
    imagesc(interB_B1_F0_1,clims)
    colorbar
    title('interB-B1-F0');
    interB_B2_F0_1=reshape(interB_B2_F0,16,8);
    interB_B2_F0_1=interB_B2_F0_1';
    clims = [aclim bclim];
    subplot(2,2,2)
    imagesc(interB_B2_F0_1,clims)
    colorbar
    title('interB-B2-F0');
    interB_B3_F0_1=reshape(interB_B3_F0,16,8);
    interB_B3_F0_1=interB_B3_F0_1';
    clims = [aclim bclim];
    subplot(2,2,3)
    imagesc(interB_B3_F0_1,clims)
    colorbar
    title('interB-B3-F0');
    intraB_F0_1=reshape(intraB_F0,16,8);
    intraB_F0_1=intraB_F0_1';
    clims = [aclim bclim];
    subplot(2,2,4)
    imagesc(intraB_F0_1,clims)
    colorbar
    title('intraB-F0');
    

    % inter-block corr for each channel for only extension mode at 45 degree
    
    % block 1 /2 /3
    for i=1:128
        interB_B1_E45(i)=corr(B1(1).E45(:,i),B1(2).E45(:,i));
        interB_B2_E45(i)=(corr(B2(1).E45(:,i),B2(2).E45(:,i))+corr(B2(1).E45(:,i),B2(3).E45(:,i))+corr(B2(2).E45(:,i),B2(3).E45(:,i)))/3;
        interB_B3_E45(i)=(corr(B3(1).E45(:,i),B3(2).E45(:,i))+corr(B3(1).E45(:,i),B3(3).E45(:,i))+corr(B3(2).E45(:,i),B3(3).E45(:,i)))/3;
    end
    
    % intra-block corr for each channel for only extension mode at 45 degree
    for i=1:128
        
        intraBAll=0;
        for j=1:2
            for jj=1:3
                intraB=(corr(B1(j).E45(:,i),B2(jj).E45(:,i))+corr(B1(j).E45(:,i),B3(jj).E45(:,i)))/2;
                intraBAll=intraBAll+intraB;
                
            end
            
        end
        
        for j=1:3
            for jj=1:3
                intraB=corr(B2(j).E45(:,i),B3(jj).E45(:,i));
                intraBAll=intraBAll+intraB;
                
            end
            
        end
        intraB_E45(i)=intraBAll/8;
    end
    
    
    figure;
    set(gcf, 'Position', [100, 100, 2400, 600]);
    suptitle(['Correlation for each channel; ',bandname{band},' band'])
    subplot(4,1,1)
    clims = [aclim bclim];
    imagesc(interB_B1_E45,clims)
    colorbar
    title('interB-B1-E45');
    
    subplot(4,1,2)
    clims = [aclim bclim];
    imagesc(interB_B2_E45,clims)
    colorbar
    title('interB-B2-E45');
    
    subplot(4,1,3)
    clims = [aclim bclim];
    imagesc(interB_B3_E45,clims)
    colorbar
    title('interB-B3-E45');
    
    subplot(4,1,4)
    clims = [aclim bclim];
    imagesc(intraB_E45,clims)
    colorbar
    title('intraB-E45');
    
    figure;
    set(gcf, 'Position', [100, 100, 2400, 1200]);
    suptitle(['Correlation for each channel; ',bandname{band},' band'])
    interB_B1_E45_1=reshape(interB_B1_E45,16,8);
    interB_B1_E45_1=interB_B1_E45_1';
    clims = [aclim bclim];
    subplot(2,2,1)
    imagesc(interB_B1_E45_1,clims)
    colorbar
    title('interB-B1-E45');
    interB_B2_E45_1=reshape(interB_B2_E45,16,8);
    interB_B2_E45_1=interB_B2_E45_1';
    clims = [aclim bclim];
    subplot(2,2,2)
    imagesc(interB_B2_E45_1,clims)
    colorbar
    title('interB-B2-E45');
    interB_B3_E45_1=reshape(interB_B3_E45,16,8);
    interB_B3_E45_1=interB_B3_E45_1';
    clims = [aclim bclim];
    subplot(2,2,3)
    imagesc(interB_B3_E45_1,clims)
    colorbar
    title('interB-B3-E45');
    intraB_E45_1=reshape(intraB_E45,16,8);
    intraB_E45_1=intraB_E45_1';
    clims = [aclim bclim];
    subplot(2,2,4)
    imagesc(intraB_E45_1,clims)
    colorbar
    title('intraB-E45');
    
    
    
    % inter-block corr for each channel for only flexion mode at 45 degree
    % block 1 /2 /3
    for i=1:128
        interB_B1_F45(i)=corr(B1(1).F45(:,i),B1(2).F45(:,i));
        interB_B2_F45(i)=(corr(B2(1).F45(:,i),B2(2).F45(:,i))+corr(B2(1).F45(:,i),B2(3).F45(:,i))+corr(B2(2).F45(:,i),B2(3).F45(:,i)))/3;
        interB_B3_F45(i)=(corr(B3(1).F45(:,i),B3(2).F45(:,i))+corr(B3(1).F45(:,i),B3(3).F45(:,i))+corr(B3(2).F45(:,i),B3(3).F45(:,i)))/3;
    end
    
    
    % intra-block corr for each channel for only flexion mode at 45 degree
    for i=1:128
        
        intraBAll=0;
        for j=1:2
            for jj=1:3
                intraB=(corr(B1(j).F45(:,i),B2(jj).F45(:,i))+corr(B1(j).F45(:,i),B3(jj).F45(:,i)))/2;
                intraBAll=intraBAll+intraB;
                
            end
            
        end
        
        for j=1:3
            for jj=1:3
                intraB=corr(B2(j).F45(:,i),B3(jj).F45(:,i));
                intraBAll=intraBAll+intraB;
                
            end
            
        end
        
        intraB_F45(i)=intraBAll/8;
    end
    
    figure;
    set(gcf, 'Position', [100, 100, 2400, 600]);
    suptitle(['Correlation for each channel; ',bandname{band},' band'])
    subplot(4,1,1)
    clims = [aclim bclim];
    imagesc(interB_B1_F45,clims)
    colorbar
    title('interB-B1-F45');
    
    subplot(4,1,2)
    clims = [aclim bclim];
    imagesc(interB_B2_F45,clims)
    colorbar
    title('interB-B2-F45');
    
    subplot(4,1,3)
    clims = [aclim bclim];
    imagesc(interB_B3_F45,clims)
    colorbar
    title('interB-B3-F45');
    
    subplot(4,1,4)
    clims = [aclim bclim];
    imagesc(intraB_F45,clims)
    colorbar
    title('intraB-F45');
    
    figure;
    set(gcf, 'Position', [100, 100, 2400, 1200]);
    suptitle(['Correlation for each channel; ',bandname{band},' band'])
    interB_B1_F45_1=reshape(interB_B1_F45,16,8);
    interB_B1_F45_1=interB_B1_F45_1';
    clims = [aclim bclim];
    subplot(2,2,1)
    imagesc(interB_B1_F45_1,clims)
    colorbar
    title('interB-B1-F45');
    interB_B2_F45_1=reshape(interB_B2_F45,16,8);
    interB_B2_F45_1=interB_B2_F45_1';
    clims = [aclim bclim];
    subplot(2,2,2)
    imagesc(interB_B2_F45_1,clims)
    colorbar
    title('interB-B2-F45');
    interB_B3_F45_1=reshape(interB_B3_F45,16,8);
    interB_B3_F45_1=interB_B3_F45_1';
    clims = [aclim bclim];
    subplot(2,2,3)
    imagesc(interB_B3_F45_1,clims)
    colorbar
    title('interB-B3-F45');
    intraB_F45_1=reshape(intraB_F45,16,8);
    intraB_F45_1=intraB_F45_1';
    clims = [aclim bclim];
    subplot(2,2,4)
    imagesc(intraB_F45_1,clims)
    colorbar
    title('intraB-F45');
    
    
    % inter-block corr for each channel for only extension mode at 90 degree
    
    % block 1 /2 /3
    for i=1:128
        interB_B1_E90(i)=(corr(B1(1).E90(:,i),B1(2).E90(:,i))+corr(B1(1).E90(:,i),B1(3).E90(:,i))+corr(B1(2).E90(:,i),B1(3).E90(:,i)))/3;
        interB_B2_E90(i)=(corr(B2(1).E90(:,i),B2(2).E90(:,i))+corr(B2(1).E90(:,i),B2(3).E90(:,i))+corr(B2(2).E90(:,i),B2(3).E90(:,i)))/3;
        interB_B3_E90(i)=(corr(B3(1).E90(:,i),B3(2).E90(:,i))+corr(B3(1).E90(:,i),B3(3).E90(:,i))+corr(B3(2).E90(:,i),B3(3).E90(:,i)))/3;
    end
    
    % intra-block corr for each channel for only extension mode at 90 degree
    for i=1:128
        
        intraBAll=0;
        for j=1:3
            for jj=1:3
                intraB=(corr(B1(j).E90(:,i),B2(jj).E90(:,i))+corr(B1(j).E90(:,i),B3(jj).E90(:,i))+corr(B2(j).E90(:,i),B3(jj).E90(:,i)))/3;
                intraBAll=intraBAll+intraB;
                
            end
            
        end
        intraB_E90(i)=intraBAll/9;
    end
    
    figure;
    set(gcf, 'Position', [100, 100, 2400, 600]);
    suptitle(['Correlation for each channel; ',bandname{band},' band'])
    subplot(4,1,1)
    clims = [aclim bclim];
    imagesc(interB_B1_E90,clims)
    colorbar
    title('interB-B1-E90');
    
    subplot(4,1,2)
    clims = [aclim bclim];
    imagesc(interB_B2_E90,clims)
    colorbar
    title('interB-B2-E90');
    
    subplot(4,1,3)
    clims = [aclim bclim];
    imagesc(interB_B3_E90,clims)
    colorbar
    title('interB-B3-E90');
    
    subplot(4,1,4)
    clims = [aclim bclim];
    imagesc(intraB_E90,clims)
    colorbar
    title('intraB-E90');
    
   
    figure;
    set(gcf, 'Position', [100, 100, 2400, 1200]);
    suptitle(['Correlation for each channel; ',bandname{band},' band'])
    interB_B1_E90_1=reshape(interB_B1_E90,16,8);
    interB_B1_E90_1=interB_B1_E90_1';
    clims = [aclim bclim];
    subplot(2,2,1)
    imagesc(interB_B1_E90_1,clims)
    colorbar
    title('interB-B1-E90');
    interB_B2_E90_1=reshape(interB_B2_E90,16,8);
    interB_B2_E90_1=interB_B2_E90_1';
    clims = [aclim bclim];
    subplot(2,2,2)
    imagesc(interB_B2_E90_1,clims)
    colorbar
    title('interB-B2-E90');
    interB_B3_E90_1=reshape(interB_B3_E90,16,8);
    interB_B3_E90_1=interB_B3_E90_1';
    clims = [aclim bclim];
    subplot(2,2,3)
    imagesc(interB_B3_E90_1,clims)
    colorbar
    title('interB-B3-E90');
    intraB_E90_1=reshape(intraB_E90,16,8);
    intraB_E90_1=intraB_E90_1';
    clims = [aclim bclim];
    subplot(2,2,4)
    imagesc(intraB_E90_1,clims)
    colorbar
    title('intraB-E90');
    
    
    
    % inter-block corr for each channel for only flexion mode at 90 degree
    % block 1 /2 /3
    for i=1:128
        interB_B1_F90(i)=(corr(B1(1).F90(:,i),B1(2).F90(:,i))+corr(B1(1).F90(:,i),B1(3).F90(:,i))+corr(B1(2).F90(:,i),B1(3).F90(:,i)))/3;
        interB_B2_F90(i)=(corr(B2(1).F90(:,i),B2(2).F90(:,i))+corr(B2(1).F90(:,i),B2(3).F90(:,i))+corr(B2(2).F90(:,i),B2(3).F90(:,i)))/3;
        interB_B3_F90(i)=(corr(B3(1).F90(:,i),B3(2).F90(:,i))+corr(B3(1).F90(:,i),B3(3).F90(:,i))+corr(B3(2).F90(:,i),B3(3).F90(:,i)))/3;
    end
    
    
    % intra-block corr for each channel for only flexion mode at 0 degree
    for i=1:128
        
        intraBAll=0;
        for j=1:3
            for jj=1:3
                intraB=(corr(B1(j).F90(:,i),B2(jj).F90(:,i))+corr(B1(j).F90(:,i),B3(jj).F90(:,i))+corr(B2(j).F90(:,i),B3(jj).F90(:,i)))/3;
                intraBAll=intraBAll+intraB;
                
            end
            
        end
        intraB_F90(i)=intraBAll/9;
    end
   
    figure;
    set(gcf, 'Position', [100, 100, 2400, 600]);
    suptitle(['Correlation for each channel; ',bandname{band},' band'])
    subplot(4,1,1)
    clims = [aclim bclim];
    imagesc(interB_B1_F90,clims)
    colorbar
    title('interB-B1-F90');
    
    subplot(4,1,2)
    clims = [aclim bclim];
    imagesc(interB_B2_F90,clims)
    colorbar
    title('interB-B2-F90');
    
    subplot(4,1,3)
    clims = [aclim bclim];
    imagesc(interB_B3_F90,clims)
    colorbar
    title('interB-B3-F90');
    
    subplot(4,1,4)
    clims = [aclim bclim];
    imagesc(intraB_F90,clims)
    colorbar
    title('intraB-F90');
    
    figure;
    set(gcf, 'Position', [100, 100, 2400, 1200]);
    suptitle(['Correlation for each channel; ',bandname{band},' band'])
    interB_B1_F90_1=reshape(interB_B1_F90,16,8);
    interB_B1_F90_1=interB_B1_F90_1';
    clims = [aclim bclim];
    subplot(2,2,1)
    imagesc(interB_B1_F90_1,clims)
    colorbar
    title('interB-B1-F90');
    interB_B2_F90_1=reshape(interB_B2_F90,16,8);
    interB_B2_F90_1=interB_B2_F90_1';
    clims = [aclim bclim];
    subplot(2,2,2)
    imagesc(interB_B2_F90_1,clims)
    colorbar
    title('interB-B2-F90');
    interB_B3_F90_1=reshape(interB_B3_F90,16,8);
    interB_B3_F90_1=interB_B3_F90_1';
    clims = [aclim bclim];
    subplot(2,2,3)
    imagesc(interB_B3_F90_1,clims)
    colorbar
    title('interB-B3-F90');
    intraB_F90_1=reshape(intraB_F90,16,8);
    intraB_F90_1=intraB_F90_1';
    clims = [aclim bclim];
    subplot(2,2,4)
    imagesc(intraB_F90_1,clims)
    colorbar
    title('intraB-F90');
    
   %----------------------------------------------------------------------- 
   % inter-block corr for each channel for extension&flexion modes at 0 degree
   % intra-block corr for each channel for extension&flexion modes at 0 degree
   
   % inter-block corr for each channel for only extension mode at 45 degree
   % intra-block corr for each channel for only extension mode at 45 degree
   
   % inter-block corr for each channel for only flexion mode at 45 degree
   % intra-block corr for each channel for only flexion mode at 45 degree
   
   % inter-block corr for each channel for extension&flexion modes at 45 degree
   % intra-block corr for each channel for extension&flexion modes at 45 degree
   
   % inter-block corr for each channel for only extension mode at 90 degree
   % intra-block corr for each channel for only extension mode at 90 degree
   
   % inter-block corr for each channel for only flexion mode at 90 degree
   % intra-block corr for each channel for only flexion mode at 90 degree
   
   % inter-block corr for each channel for extension&flexion modes at 90 degree
   % intra-block corr for each channel for extension&flexion modes at 90 degree
   
   
   % inter-block corr for each channel for only extension mode for all degree
   % intra-block corr for each channel for only extension mode for all degree
   
   % inter-block corr for each channel for only flexion mode for all degree
   % intra-block corr for each channel for only flexion mode for all degree
   
   % inter-block corr for each channel for extension&flexion modes for all degree
   % intra-block corr for each channel for extension&flexion modes for all degree
end
 %%   
 %%%%%  % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 


