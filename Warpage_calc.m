%% User inputs
sideImagePath = '';      % if empty use synthetic example
px_to_mm = 0.1;          % calibration from pixel to mm 
polyDeg = 2;             % for quadratic  

%% Load or generate side image
if ~isempty(sideImagePath) && isfile(sideImagePath)
    sideI = imread(sideImagePath);
else
    W=900; H=300; x=linspace(-1,1,W); k=0.8;
    y = 140 - (k*(x.^2)*60);
    mask = false(H,W);
    for c=1:W
        t = max(1,min(H,round(y(c))));
        mask(t:H-30,c)=true;
    end
    sideI = uint8(~mask)*255;
end

%% Silhouette extraction
if size(sideI,3)==3, Igray=rgb2gray(sideI); else, Igray=sideI; end %grayscale conversion
Igray = im2double(Igray);
bw = imfill(bwareaopen(Igray < graythresh(Igray),300),'holes');
[Himg,Wimg] = size(bw);

%% Top-edge sampling (extracting edges)
topY = nan(1,Wimg);
for c=1:Wimg
    ys = find(bw(:,c));
    if ~isempty(ys), topY(c)=min(ys); end
end
valid = ~isnan(topY);
x_pix = find(valid);
y_pix = topY(valid);

%% Fit & convert
p = polyfit(x_pix,y_pix,polyDeg);
y_fit = polyval(p,x_pix);
nEdge = max(1,round(numel(y_pix)*0.05));
baseline = median([median(y_pix(1:nEdge)), median(y_pix(end-nEdge+1:end))]);
z_mm = (baseline - y_fit) * px_to_mm;

%% Metrics
maxLift = max(z_mm);
pptp = max(z_mm)-min(z_mm);
rmsWarp = sqrt(mean((z_mm-mean(z_mm)).^2));

%% Curvature estimate (for quadratic)
if polyDeg==2
    x_mm = (x_pix-mean(x_pix))*px_to_mm;
    p_mm = polyfit(x_mm,z_mm,2);
    a = p_mm(1);
    k_m = 2*a*1000;
    if k_m~=0
        R = 1/k_m;
        L = (max(x_mm)-min(x_mm))/1000;
        delta = 1000*(R*(1-cos(L/(2*R))));  %wang's model prediction
    else
        R = Inf; delta = NaN;
    end
else
    k_m = NaN; R=NaN; delta=NaN;
end

%% Plots
figure; imshow(sideI); hold on;
plot(x_pix,y_pix,'.g'); plot(x_pix,y_fit,'r','LineWidth',2);
title('Silhouette fit');

figure; plot((x_pix-mean(x_pix))*px_to_mm, z_mm,'.-b');
xlabel('x (mm)'); ylabel('z (mm)'); grid on;

%% Output
fprintf('\nMax lift: %.3f mm\n', maxLift);
fprintf('Peak-to-peak: %.3f mm\n', pptp);
fprintf('RMS warp: %.3f mm\n', rmsWarp);
fprintf('Curvature: %.4e 1/m,  R=%.3f m,  delta≈%.3f mm\n', k_m, R, delta);
