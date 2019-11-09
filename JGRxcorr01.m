clear all; clc; warning off; close all;

%% add path
addpath(genpath('C:\Users\DELL\OneDrive - OBA GG\11_yang_code_Linux\05_matlab\SeismicLab\codes'));
addpath(genpath('C:\Users\DELL\OneDrive - OBA GG\09_Linux_Software\Software\TimeFreq\Wavelet\YangWavelet'));
addpath matdata

flag = 'nnyn';

load data_mario_proc

% data_demo = data(:,end/2-800:end/2+800)';

data_demo = data;

offset = linspace(0,size(data_demo,2)*1,size(data_demo,2));
depth = offset;
xcdp = linspace(5,50,length(depth));

csg1 = data_demo;
vini = 2000/730;
%% define analysis parameters

wlen = 64;                        % window length (recomended to be power of 2)
dt = 0.004;
dj = 0.05; % resoltuions
s0 = 0.001;
%s0 = 0.0118;
ub = 100;
lb = 1;
Nr = 7;
Nc = 3;
nt = size(csg1,1);
ntrc = size(csg1,2);
nti = nt*2-1;
nfft = nti;
win = 25;
rec = size(csg1,2);
FilterT = 12;
FilterF = 5;
j1=ceil((log(length(csg1(:,1))*dt/s0)/log(2))/dj);
time = linspace(0,length(data_demo)*dt,length(data_demo));

%%
for iRec = 1:rec
    %for iRec = 1:1
    
    up = csg1(:,iRec);
    
    %% compute the CWT
    % [waveup, periodup, scaleup, coiup, djup,paramoutup, kup] = contwt(up,dt,[],dj,s0,[],'MORLET',6);
    [waveup, periodup, scaleup, coiup, djup,paramoutup, kup] = contwt(up,dt,[],dj,s0);
    
    
    % do wavelet domain 2D filteringfi
    % wavex2  = wavex(:,nti:end);
    wavex2  = waveup;
    amp = abs(wavex2);
    pha = angle(wavex2);
    amp4 = amp;
    
    %% output
    %amp4 = amp;
    
    
    [framefil,label]= ImgSeg02(amp4,iRec,vini,distance);
    wavex1 = framefil.*complex(cos(pha),sin(pha));
    
    %% apply wavelet domain filter and sum over desired frequency
    wavexrec(:,iRec) = invcwt(wavex1, 'MORLET', scaleup, paramoutup, kup);
    wavexCuve(:,:,iRec) = wavex1;
    
    %% save the frames
    labelCubeAll(:,:,iRec) = label;
    labelCubePick(:,:,iRec) = framefil;
    Xcube(:,:,iRec) = amp;
    
%     figure(11);
%     subplot(4,1,[1 3]);
%     imagesc(amp);
%     subplot(4,1,4);
%     plot(csg1(:,iRec));
%     drawnow; pause(2);
end

%% normalization
for itrc = 1:size(csg1,2)
    data_demo(:,itrc) = data_demo(:,itrc)./max(abs(data_demo(:,itrc)));
    wavexrec(:,itrc) = wavexrec(:,itrc)./max(abs(wavexrec(:,itrc)));
end


%%
close all;
% create ticks and reda
timeticks = linspace(1, size(data_demo, 2), numel(time));
dt = 0.0008;
nt = size(Xcube,2);
time = linspace(nt*dt/2*(-1),nt*dt/2,nt)*1000;

figure(30);
subplot(2,1,1);
wigb(data_demo,0.7,offset,time);
ylabel('Time (ms)');
title('(a)', 'FontSize', 15);
hold on;
[maxlin,idx1]=max((data_demo));
%plot(offset, time(idx1),'r','LineWidth',2);
hold off;
%set(gca, 'YTick', timeticks(1:50:end), 'YTickLabel', round(time(1:50:end)));


subplot(2,1,2);
wigb(wavexrec,0.7,offset,time);
title('(b)', 'FontSize', 15);
hold on;
[maxlin,idx2]=max(diff(wavexrec));
%plot(offset,smooth(time(idx2)),'r','LineWidth',2);
hold off;
xlabel('Station Number');
ylabel('Time (ms)');
%set(gca, 'YTick', timeticks(1:50:end), 'YTickLabel', round(time(1:50:end)));


%save('LabelMatrix02.mat','labelCubeAll');
%save('OutputData03.mat','labelCubePick');
%save('WXcorrData01.mat','data_demo','Xcube','wavexrec');
before_process = data_demo;
after_process = wavexrec;



%% draw professional plot
figure(31);
subplot(1,2,1);
    
for iRec = 1:rec
    plot(before_process(:,iRec)*20+distance(iRec),'k');    
    hold on
    drawnow;
end

xlim([0 800])
set(gca,'XDir','normal','YDir','reverse');
xlabel('Time (s)');
ylabel('Distante (km)')
title('(a)', 'FontSize', 15);


subplot(1,2,2);    
for iRec = 1:rec
    plot(after_process(:,iRec)*20+distance(iRec),'k');    
    hold on
    drawnow;
end

xlim([0 800])
set(gca,'XDir','normal','YDir','reverse');
xlabel('Time (s)');
ylabel('Distante (km)')
title('(b)', 'FontSize', 15);
saveFigure(31,'Crosscorrelation vs WaveletKalman',flag);
save('data_yang_zhao','before_process','after_process','wavexCuve');