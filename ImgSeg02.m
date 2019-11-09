function [framefil L]= ImgSeg02(frame,rec,vini,distance)
%% image segmentation
% step 1: origianl data downsampled in time domain
fig = 1;
sampleN = 100;
I = frame;
thresH = 0.15;
threshold = thresH*max(I(:));
I(I <= threshold) = 0;
I = exp(bsxfun(@rdivide,I,max(I(:))));

% step 2: mark the foreground object
% Opening-by-reconstruction
se = strel('disk', 10);
Ie = imerode(I, se);
Iobr = imreconstruct(Ie, I);

figure(14);
% imagesc(time,freq,flipud(Iobr));
imagesc((Iobr));
set(gca,'YDir','normal');
hcb=colorbar;
title(hcb,'Magnitude');
xlabel('Time (s)');
ylabel('Frequency  (Hz) ');

% Closing-by-reconstruction
Iobrd = imdilate(Iobr, se);
Iobrcbr = imreconstruct(imcomplement(Iobrd), imcomplement(Iobr));
Iobrcbr = imcomplement(Iobrcbr);
%     figure(4); clf;
%     imagesc(Iobrcbr);colorbar;
%     title('Opening-closing by reconstruction (Iobrcbr)');

fgm = imregionalmax(Iobrcbr);
% figure(100);
% subplot(2,3,2);
% imagesc(time,freq,flipud(fgm));
% title('Regional maxima of opening-closing by reconstruction');
% set(gca,'YDir','normal');


% label area and choose low frequency label
L = bwlabel(fgm);

% predicted wave-mode location based on wave velocity, vini, rec #
predf = 125;
predt1 = distance(rec)/vini;
%predt2 = size(fgm,2)/2 + distance/vini;
%predt1
% find two peaks if needed
display(distance(rec))

for ilabel = 1:4
    [r, c] = find(L==ilabel);
    dist1(ilabel) = sqrt(abs(median(c)-predt1).^2+abs(median(r)-predf).^2);
    rmed1(ilabel) = mean(r);
    
    %dist2(ilabel) = sqrt(abs(median(c)-predt2).^2+abs(median(r)-predf).^2);    
end

% find the label with lowest frequency
[dum,labelN1] = min(dist1);
%[dum,labelN2] = min(dist2);

% extract ffrom lowest frequency label and construct mask
Lren = L; %(1:sampleN:end,:);
Lds1 = zeros(size(Lren));
Lds2 = zeros(size(Lren));

framemask1 = frame;
framemask1(find(L~=labelN1))=0;
Lds1(find(Lren==labelN1))=1;
Ldsmth1 = smooth2a(Lds1,10,10);
Ldsmth1 = Ldsmth1/max(max(Ldsmth1));

% framemask2 = frame;
% framemask2(find(L~=labelN2))=0;
% Lds2(find(Lren==labelN2))=1;
% Ldsmth2 = smooth2a(Lds2,10,10);
% Ldsmth2 = Ldsmth2/max(max(Ldsmth2));

% Ldsmth = Ldsmth2 + Ldsmth1;
Ldsmth = Ldsmth1;
framefil = frame.*Ldsmth;

%%

figure(12);
imagesc((L));
%title('Label matrix');
set(gca,'YDir','normal');
hcb=colorbar;
title(hcb,'Magnitude');
xlabel('Time (s)');
ylabel('Frequency  (Hz) ');

hold on;
plot(predt1,predf,'y*');
%hold on;
%plot(predt2,predf,'yo');
hold off;

figure(13);
subplot(4,1,[1 3]);
imagesc((framefil));
title('After Mask');
set(gca,'YDir','normal');
hcb=colorbar;
title(hcb,'Magnitude');
ylabel('Frequency  (Hz) ');
title('(a)', 'FontSize', 15);

return
