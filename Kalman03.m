clear all; clc; warning off; close all;

%% add path
addpath(genpath('C:\Users\zhaoy\OneDrive - OBA GG\11_yang_code_Linux\05_matlab\SeismicLab\codes'));
addpath(genpath('C:\Users\zhaoy\OneDrive - OBA GG\09_Linux_Software\Software\TimeFreq\Wavelet\YangWavelet'));
flagfig = 'nnyn';
addpath data

%% load Wavelet Xcorr Data Cube and label matrix ( just for presentation illustration)
load LabelMatrix02.mat
load OutputData03.mat
load WXcorrData01.mat

dt = 0.0008;
dx = 10;
dz = 10;

%wavexrec = wavexrec(2000:4000,:);
xcsg = data_demo;
depth = linspace(0,1195*10,1195);
offset = linspace(1990,1990+size(xcsg,2)*20,size(xcsg,2));
depth = offset;

nt = size(Xcube,2);
time = linspace(nt*dt/2*(-1),nt*dt/2,nt)*1000;
freq = linspace(1,40,size(Xcube,1));

%% Kalman filter initialization
R=[[0.2845,0.0045]',[0.0045,0.0455]'];
% H=[[1,0]',[0,1]',[0,0]',[0,0]'];
H = eye(4);
Q=0.01*eye(4);
P = 100*eye(4);
dt=1;
A=[[1,0,0,0]',[0,1,0,0]',[dt,0,1,0]',[0,dt,0,1]'];
g = 6; % pixels^2/time step
Bu = [0,0,0,g]';
kfinit=0;
x=zeros(10,4);
MC = size(labelCubeAll,1);
MR = size(labelCubeAll,2);
scrsz = get(0,'ScreenSize');

%% loop over all image frames
for iframe = 1:1:size(labelCubeAll,3)
    %for iframe = 1:2000
    %for iframe = [30 160 400 800]
    % load image
    Imwork = labelCubePick(:,:,iframe);
    ImPlot = labelCubeAll(:,:,iframe);
    
    
    % extract the center and box edges
    stats = regionprops(Imwork,['basic']);
    centroid = stats(1).Centroid;
    radius = sqrt(stats(1).Area/pi);
    BoundingBox = stats(1).BoundingBox;
    Bounding(iframe,:) = BoundingBox;
    cc(iframe) = centroid(1);
    cr(iframe) = centroid(2);
    flag = 1;
    
    %[cc(iframe),cr(iframe),radius,flag] = extractball(Imwork,Imback,iframe);
    if flag==0
        continue
    end
    
    %     hold on
    %     for c = -1*radius: radius/20 : 1*radius
    %         r = sqrt(radius^2-c^2);
    %         plot(cc(iframe)+c,cr(iframe)+r,'g.')
    %         plot(cc(iframe)+c,cr(iframe)-r,'g.')
    %     end
    % Kalman update
    iframe
    if kfinit==0
        xp = [MC/2,MR/2,0,0]';
    else
        xp=A*x(iframe-1,:)' + Bu;
    end
    kfinit=1;
    PP = A*P*A' + Q;
    K = PP*H'*inv(H*PP*H');
    % K = PP*H'*inv(H*PP*H'+R);
    % x(iframe,:) = (xp + K*([cc(iframe),cr(iframe)]' - H*xp))';
    x(iframe,:) = (xp + K*(BoundingBox' - H*xp))';
    P = (eye(4)-K*H)*PP;
    
    hold off;
    figure(66);clf;
    %set(gcf,'Position',[1 scrsz(4)*0.9 scrsz(3)*0.9 scrsz(4)*0.9]);
    subplot(2,2,1);
    wigb(xcsg,0.5,offset,time);  
    hold on;
    line([offset(iframe) offset(iframe)],[-1 1],'Color','red','LineStyle','-','LineWidth',12);
    hold off;
    xlabel('Measured depth (m)');
    title('(a)', 'FontSize', 15);
    ylabel('Time (ms)');

    ax2=subplot(2,2,2);
    imagesc(time,freq,flipud(abs(Xcube(:,:,iframe))));
    ylabel('Frequency (Hz) ');
    set(gca,'YDir','normal');
    xlabel('Time (ms)');
    colormap('default');
    title('(b)', 'FontSize', 15);
    
    subplot(2,2,4);
    imagesc(ImPlot);
    hold on;
    h = rectangle('Position',x(iframe,:),'LineStyle',':','LineWidth',4);
    set(h,'EdgeColor','red');
    colormap('default');
    ylabel('Frequency (Hz) ');
    %set(gca,'YDir','normal');
    xlabel('Time (ms)');
    title('(c)', 'FontSize', 15);
    
    % create ticks and reda
    timeticks = linspace(1, size(ImPlot, 2), numel(time));
    set(gca, 'XTick', timeticks(1:50:end), 'XTickLabel', round(time(1:50:end)));
    
    freqticks = linspace(1, size(ImPlot, 1), numel(fliplr(freq)));
    set(gca, 'YTick', freqticks(1:25:end), 'YTickLabel', round(fliplr(freq(1:25:end))));   
    hold off;
    
    subplot(2,2,3);
    tmp = flipud(abs(Imwork));
    frametmp = tmp/max(max(abs(tmp)))*2.4;
    
    imagesc(time,freq,frametmp);
    ylabel('Frequency (Hz) ');
    set(gca,'YDir','normal');
    xlabel('Time (ms)');
    colormap('default');
    title('(d)', 'FontSize', 15);
    
    hcb=colorbar;
    title(hcb,'Magnitude');
    set(hcb,'position',[0.95 0.25 0.01 0.5]);
    
    %if isequal(iframe,30) ||  isequal(iframe,130) ||  isequal(iframe,160) ||  isequal(iframe,400) || isequal(iframe,800)
        saveName = strcat('tracking_kalman_frame',num2str(iframe));
        saveFigure(66,saveName,flagfig);
    %end
    
    drawnow;
    set(gcf,'color','w');
    %KalmanMov(iframe)=getframe(gcf);
    pause(1)
end

%%
figure(62);
subplot(2,1,1);
plot(offset,(x(:,1)-cc')/20,'b');
title('(a)', 'FontSize', 15);
grid on;
ylabel('Frequency (Hz) ');
    

subplot(2,1,2)
plot(offset,(x(:,2)-cr')*0.0008,'b');
xlabel('Measured depth (m)');
ylabel('Time (ms)');
title('(b)', 'FontSize', 15);
grid on;
%legend('predict time','predict freq','measure time','measure freq');
saveFigure(62,'Predict_vs_Measure_microseismic',flagfig);