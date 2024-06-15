%% Load DICOM files

clear all
close all
clc

%Choose path and get filenames
filePath = uigetdir();
cd(filePath)
dicomData = dir('**/*.dcm'); %list of all dicom files in the folder

dicomNames = extractfield(dicomData,'name'); %now all dicomNames in a cell
dicomNames = string(dicomNames); %convert it into strings
numberOfFiles = length(dicomNames);

%Load all dicomfiles
dicom1 = dicomread(dicomNames(1));
CTA = zeros(size(dicom1,1),size(dicom1,2),numberOfFiles); %predefine CTA matrix
dicomNames_sort = sort(dicomNames,'descend'); %sort dicomfiles in descending order
for i = 1:numberOfFiles
    CTA(:,:,i) = dicomread(dicomNames_sort(i));
end

%Hounsfield units terugzetten: 
dicom_information = dicominfo(dicomNames(1));
CTA = dicom_information.RescaleIntercept + dicom_information.RescaleSlope * CTA; 
originalCTA = CTA;
cd('C:\Users\Astrid\surfdrive\Documents\MATLAB\Carotis segmentatie'); %go to correct folder
addpath(genpath(cd)); %add current directory to path

%% Crop CTA dataset to region of interest
% Not every CTA consists of an equal number of slices. Take as lower bound
% the slice on 45% and as upper bound the slice on 70%. 
lowerbound = round(0.45*size(originalCTA,3));
upperbound = round(0.7*size(originalCTA,3));
CTA = CTA(:,:,lowerbound:upperbound);
CTA = flip(CTA,3); %flip 3rd dimension of CTA such that slice 1 is the lowest slice (containing CCA)

%% Make binary of whole CTA dataset: ---------------------------------------------THRESHOLDING-----------------------------------------
%lumenMin_temp = 220;         % Define thresholds 
%lumenMax_temp = 445;         % Define thresholds 
% calciumMin = 426;       % Define thresholds 
% calciumMax = 1224;      % Define thresholds

%binaryCTA_temp = (CTA >= lumenMin_temp) & (CTA <= lumenMax_temp); % Create mask

% Select manually carotid in first slice
figure('Position', [100 300 600 600]);
imshow(CTA(:,:,1),[]); %show orinal 1st slice
annotation('textbox', [0.08 0.5 0.8 0.1], ...
    'String', 'R', ...
    'Color', [0 0 0], ...
    'FontWeight', 'bold', ...
    'EdgeColor', 'none')
annotation('textbox', [0.92 0.5 0.8 0.1], ...
    'String', 'L', ...
    'Color', [0 0 0], ...
    'FontWeight', 'bold', ...
    'EdgeColor', 'none')
% figure('Position', [800 300 600 600]);
% imshow(binaryCTA_temp(:,:,1),[]); %show masked 1st slice
% annotation('textbox', [0.08 0.5 0.8 0.1], ...
%     'String', 'R', ...
%     'Color', [0 0 0], ...
%     'FontWeight', 'bold', ...
%     'EdgeColor', 'none')
% annotation('textbox', [0.92 0.5 0.8 0.1], ...
%     'String', 'L', ...
%     'Color', [0 0 0], ...
%     'FontWeight', 'bold', ...
%     'EdgeColor', 'none')
title('Click on the carotid artery in this image and hit enter')
[x,y] = getpts; %click on the vessel and hit enter
close all %close images

x_CCA = x; %x coordinate common carotid artery
y_CCA = y; %y coordinate common carotid artery

app_running = app_thresholdselection;%launch app for threshold selection
waitfor(app_running)  % wait for app to close

binaryCTA = (CTA >= lumenMin_app) & (CTA <= lumenMax_app); % Create mask


%% Keep the connected components
 
CC = bwconncomp(binaryCTA); %find all connected components in binaryCTA
labeledCC = labelmatrix(CC); %er vanuitgaande dat de carotis in elke slice tot hetzelfde 
carotidLabel = labeledCC(ceil(y_CCA(1)),ceil(x_CCA(1)),1); %index of label where carotid is located
%carotidIdx = find(labeledCC==carotidLabel); %index of pixels that belong to carotid
[R,C,P] = ind2sub(size(CTA),find(labeledCC==carotidLabel)); %row, column and plane coordinates of pixels that belong tot carotid
%Now: make binary of CTA where the pixels R, C and P are 1 and rest is zero
zerosCTA = zeros(size(CTA)); %Make zero matrix of size CTA
for i=1:size(R)
    zerosCTA(R(i),C(i),P(i))=1; % Make all pixels of R, C and P ones
end
binaryCarotid = zerosCTA;

%% Smooth binary carotid (per slice)
for i=1:size(binaryCarotid,3)
    currentBinary = binaryCarotid(:,:,i);
    
    % Take boundaries for current binary image:
    [B,L] = bwboundaries(currentBinary);
    boundaries{i} = B; %save all boundaries in matrix
    
    % Blur image and threshold at 0.5
    width = 3;
    kernel = ones(width) / width^2;
    blurryImage = conv2(double(currentBinary), kernel, 'same');
    smoothImage = blurryImage > 0.5;
    
    % Take boundaries for current smoothed binary image:
    [B_smooth,L_smooth,N_smooth] = bwboundaries(smoothImage);
    boundaries_smooth{i} = B_smooth; %save all boundaries in matrix
    
    % Put smoothImage in matrix:
    binaryCarotid_smooth(:,:,i) = smoothImage;
end
    
%% Make centerlines
% Make two centerlines: CCA-ICA and CCA-ECA

% All centerpoints:
for i = 1:size(binaryCarotid_smooth,3)
    s = regionprops(binaryCarotid_smooth(:,:,i),'centroid','area','EquivDiameter');
    centerpoint = cat(1,s.Centroid);
    area = cat(1,s.Area);
    equivDiam = cat(1,s.EquivDiameter);
    centerpoints{i} = centerpoint;
    areas{i} = area;
    equivDiams{i} = equivDiam;
end

app2_running = app_imageslide4; %launch app
waitfor(app2_running) % wait for app to close

% Select ICA at latest slice:
% figure, imshow(binaryCarotid_smooth(:,:,size(binaryCarotid_smooth,3)));
% title('Click on the internal carotid artery in this image and hit enter')
% [x_ICA,y_ICA] = getpts; %click on the vessel and hit enter
% close all %close image
% 
% Select ECA at latest slice:
% figure, imshow(binaryCarotid_smooth(:,:,size(binaryCarotid_smooth,3)));
% title('Click on the external carotid artery in this image and hit enter')
% [x_ECA,y_ECA] = getpts; %click on the vessel and hit enter
% close all %close image

%--------------------------------------------------------------------------
%HIER IETS INBOUWEN VOOR NO ICA EN/OF NO ECA
%--------------------------------------------------------------------------

%App delivers point_ICA and point_ECA
x_ICA = point_ICA(:,1);
y_ICA = point_ICA(:,2);
x_ECA = point_ECA(:,1);
y_ECA = point_ECA(:,2);
%App also delivers sliceNr_ICA and sliceNr_ECA. This is the new amount of
%slices

%App also delivers sliceNr_sten. This is an string containing one or more
%slice numbers where the stenosis/stenoses is/are located. The strings need
%to be separated and converted to numeric values.
sliceNr_sten = strsplit(sliceNr_stenApp,';');
for i=1:length(sliceNr_sten)
    sten(i) = str2num(sliceNr_sten{i});
end

% %Launch app for bifurcation selection:
% app3_running = app_selectBifurcation2;
% waitfor(app3_running); %wait for app to close
% 
% %App delivers point_bifurcation 
% slice_bif = upperbound - fix(point_bifurcation(2)); %change slice nr from originalCTA to CTA coordinates

%aantal slices (n_slices) is gelijk laagste slicenummer waar ICA of ECA
%gekozen is. 
% if isempty(sliceNr_ICA) && isempty(sliceNr_ECA)
%     n_slices = 50; %--------------------------------------------------------DIT NOG INVULLEN!!!
% elseif isempty(sliceNr_ICA)
%     n_slices = sliceNr_ECA;
% elseif isempty(sliceNr_ECA)
%     n_slices = sliceNr_ICA;
% elseif sliceNr_ICA > sliceNr_ECA
%     n_slices = sliceNr_ECA;
% elseif sliceNr_ECA > sliceNr_ICA
%     n_slices = sliceNr_ICA; 
% elseif sliceNr_ICA == sliceNr_ECA
%     n_slices = sliceNr_ICA;
% end

n_slices_ICA = sliceNr_ICA;
n_slices_ECA = sliceNr_ECA;

clear center_ICA
center_ICA(n_slices_ICA,:) = [x_ICA,y_ICA]; %take manual input as 'first' centerpoint

for i=1:n_slices_ICA-1
    sliceNr = n_slices_ICA+1-i;
    centerpoints_next = centerpoints{n_slices_ICA-i}; %coordinates of centerpoint on next slice
    
%     % If there is no segmentation (occlusion) on next slice:
%     if isempty(centerpoints_next) %if there is no centerpoint (occlusion)
%         centerpoints_next = centerpoints{n_slices_ICA}; %take current centerpoint
%         occlusion_ICA = true;
%     end
    
    % Calculate distances between centerpoints
    distances = pdist2(center_ICA(n_slices_ICA+1-i,:),centerpoints_next,'euclidean'); %calculate distances from centerpoint to all centerpoints on next slice
    [mindistance,mindis_idx] = min(distances); %smallest distance & index of smallest distance
    center_ICA(n_slices_ICA-i,:) = centerpoints_next(mindis_idx,:); %closest center on next slice
    
%     equivDiam = equivDiams{sliceNr}(mindis_idx);
%     % If there is no centerpoint within distance <0.5 of the equivalent diameter:
%     if mindistance < 0.5*equivDiams{n_slices_ICA-1}
%         occlusion_ICA = true;
%         center_ICA(n_slices_ICA-1,:) = centerpoints{n_slices_ICA}; %take current centerpoint
%     else
%         center_ICA(n_slices_ICA-i,:) = centerpoints_next(mindis_idx,:); %closest center on next slice
%     end
       
end
center_ICA = [center_ICA (1:n_slices_ICA)']; %add slice numbers to centerline coordinates

% Boundaries ICA
for i=1:n_slices_ICA
    current_centerpoint = centerpoints{i};
    if size(current_centerpoint,1)>1
        for j=1:size(current_centerpoint,1)
            distance_cp(j,:) = pdist2(current_centerpoint(j,:),center_ICA(i,1:2),'euclidean');
        end
        [mindistance,mindis_idx] = min(distance_cp); %smallest distance & index of smallest distance
        boundary_ICA = boundaries_smooth{i}{mindis_idx};
        boundaries_ICA{i} = boundary_ICA;
        clear distance_cp
    elseif isempty(boundaries_smooth{i}) %in case of occlussion
        boundaries_ICA{i} = boundaries_smooth{i}; %boundaries_ICA{i} is empty
    else
        boundaries_ICA{i} = boundaries_smooth{i}{1};
    end
    clear current_centerpoint
end

clear center_ECA
center_ECA(n_slices_ECA,:) = [x_ECA,y_ECA]; %take manual input as 'first' centerpoint
for i=1:n_slices_ECA-1
    sliceNr = n_slices_ECA+1-i;
    centerpoints_next = centerpoints{n_slices_ECA-i}; %coordinates of centerpoint on next slice
    
    % If there is no segmentation (occlusion) on next slice:
    if isempty(centerpoints_next) %if there is no centerpoint (occlusion)
        centerpoints_next = centerpoints{n_slices_ECA}; %take current centerpoint
        occlusion_ECA = true;
    end
    
    distances = pdist2(center_ECA(n_slices_ECA+1-i,:),centerpoints_next,'euclidean'); %calculate distances from centerpoint to all centerpoints on next slice
    [mindistance,mindis_idx] = min(distances); %smallest distance & index of smallest distance
    center_ECA(n_slices_ECA-i,:) = centerpoints_next(mindis_idx,:); %closest center on next slice
end
center_ECA = [center_ECA (1:n_slices_ECA)']; %add slice numbers to centerline coordinates

% Boundaries ECA
for i=1:n_slices_ECA
    current_centerpoint = centerpoints{i};
    if size(current_centerpoint,1)>1
        for j=1:size(current_centerpoint,1)
            distance_cp(j,:) = pdist2(current_centerpoint(j,:),center_ECA(i,1:2),'euclidean');
        end
        [mindistance,mindis_idx] = min(distance_cp); %smallest distance & index of smallest distance
        boundary_ECA = boundaries_smooth{i}{mindis_idx};
        boundaries_ECA{i} = boundary_ECA;
        clear distance_cp
    elseif isempty(boundaries_smooth{i}) %in case of occlussion
        boundaries_ECA{i} = boundaries_smooth{i}; %boundaries_ECA{i} is empty
    else boundaries_ECA{i} = boundaries_smooth{i}{1};
    end
    clear current_centerpoint
end

%% Calculate diameters ICA and other properties on orthogonal slices

for i=1:length(boundaries_ICA)
    % select correct vessel on orthogonal slice. get binary of only ICA:
    orthSlice = binaryCarotid_smooth(:,:,i); %select current orthogonal slice
    [B,L] = bwboundaries(orthSlice); %find boundaries and labels of all regions
    % Find the label number where centerpoint is located:
    L_ICA = L(round(center_ICA(i,2)),round(center_ICA(i,1)));
    % If label == 0, then find nearest white pixel and give
    % binaryCarotid_smooth that label:
    if L_ICA==0
        %calculate distance to all white pixels and choose pixels with
        %smallest distance        
        [x,y] = find(orthSlice==1); % find all nonzeropixels
        nonzeropixels = [x,y]; % array with non-zero pixels
        distances = sqrt((x - center_ICA(i,2)).^2 + ...
                            (y - center_ICA(i,1)).^2); % distances to nonzero pixels:
        [~,sdidx] = min(distances); %linear index of pixel with smallest distance 
        sd = nonzeropixels(sdidx,:); %row-column index of pixel with smallest distance
        L_ICA2 = L(sd(1),sd(2)); % label value of the pixel with smallest distance
        binary_ICA = L==L_ICA2; %binary with only the label where pixel with smallest distance is located 
    else  
        % Get binary with only the label where centerpoint is located:
        binary_ICA = L==L_ICA;
    end
    binaries_ICA{i} = binary_ICA; %save binary orthogonal images of ICA in structure
    
    boundary_orth_ICA = boundaries_ICA{i}; %select current boundary
    centroid = [center_ICA(i,2) center_ICA(i,1)]; %coordinates of current centroid
    
    %Show binary with boundary and center:
%     figure, imshow(binary_ICA);
%     hold on;
%     plot(center_ICA(i,1),center_ICA(i,2),'r.');
%     plot(boundary_orth_ICA(:,2),boundary_orth_ICA(:,1),'b-');
    
    % For each boundary point find distance to boundary point at opposite
    % side and choose closest boundary point for calculation of diameter:
    clear diameter    
    for j=1:(0.5*length(boundary_orth_ICA))
        boundaryj = boundary_orth_ICA(j,:); %current boundarypoint

        %search for closest point on boundary:
        index = j+floor(0.25*length(boundary_orth_ICA)):(j+floor(0.25*length(boundary_orth_ICA))+floor(0.5*length(boundary_orth_ICA))); %indexes of search area
        boundary_enlarged = repmat(boundary_orth_ICA,2,1); %place boundaries matrix after boundaries matrix for making 'a loop' while choosing the right diameter
        boundary_search = boundary_enlarged(index,:); %subset of boundary matrix, only values relevant for search
        distances = point_to_line_distance(boundary_search, boundaryj, centroid); %distances from line boundaryj-centroid to the points of boundary_search 

        [val,ind] = min(distances); %find value and index of boundarypoint closest to line 
        diameter(j) = pdist([boundaryj;boundary_search(ind,:)]); %Calculate diameter using eucledian distance:
        diameter_plot_ICA{i}(j,:) = [boundaryj boundary_search(ind,:)]; %Save coordinates of diameter
    end
    
    %regionprops:
    stats_ICA = regionprops(binary_ICA,'Area','Circularity','Eccentricity','EquivDiameter','MajorAxisLength','MinorAxisLength','Orientation','Centroid');
    
    diameters_ICA{i} = diameter; %save all diameter values per slice
    mindiameters_ICA(i) = min(diameters_ICA{i}); %smallest diameter for each slice
    maxdiameters_ICA(i) = max(diameters_ICA{i}); %largest diameter for each slice
    avgdiameters_ICA(i) = mean(diameters_ICA{i}); %average diameter for each slice
    area_ICA(i) = pi*((0.5*avgdiameters_ICA(i))^2); %area of vessel for each slice
         
    s_ICA(i)= stats_ICA;
end


%% Calculate diameters ECA and other properties on orthogonal slices

for i=1:length(boundaries_ECA)
    % select correct vessel on orthogonal slice. get binary of only ECA:
    orthSlice = binaryCarotid_smooth(:,:,i); %select current orthogonal slice
    [B,L] = bwboundaries(orthSlice); %find boundaries and labels of all regions
    % Find the label number where centerpoint is located:
    L_ECA = L(round(center_ECA(i,2)),round(center_ECA(i,1)));
    % If label == 0, then find nearest white pixel and give
    % binaryCarotid_smooth that label:
    if L_ECA==0
        %calculate distance to all white pixels and choose pixels with
        %smallest distance        
        [x,y] = find(orthSlice==1); % find all nonzeropixels
        nonzeropixels = [x,y]; % array with non-zero pixels
        distances = sqrt((x - center_ECA(i,2)).^2 + ...
                            (y - center_ECA(i,1)).^2); % distances to nonzero pixels:
        [~,sdidx] = min(distances); %linear index of pixel with smallest distance 
        sd = nonzeropixels(sdidx,:); %row-column index of pixel with smallest distance
        L_ECA2 = L(sd(1),sd(2)); % label value of the pixel with smallest distance
        binary_ECA = L==L_ECA2; %binary with only the label where pixel with smallest distance is located 
    else  
        % Get binary with only the label where centerpoint is located:
        binary_ECA = L==L_ECA;
    end
    binaries_ECA{i} = binary_ECA; %save binary orthogonal images of ECA in structure
    
    boundary_orth_ECA = boundaries_ECA{i}; %select current boundary
    centroid = [center_ECA(i,2) center_ECA(i,1)]; %coordinates of current centroid
    
    %Show binary with boundary and center:
%     figure, imshow(binary_ICA);
%     hold on;
%     plot(center_ICA(i,1),center_ICA(i,2),'r.');
%     plot(boundary_orth_ICA(:,2),boundary_orth_ICA(:,1),'b-');
    
    % For each boundary point find distance to boundary point at opposite
    % side and choose closest boundary point for calculation of diameter:
    clear diameter    
    for j=1:(0.5*length(boundary_orth_ECA))
        boundaryj = boundary_orth_ECA(j,:); %current boundarypoint

        %search for closest point on boundary:
        index = j+floor(0.25*length(boundary_orth_ECA)):(j+floor(0.25*length(boundary_orth_ECA))+floor(0.5*length(boundary_orth_ECA))); %indexes of search area
        boundary_enlarged = repmat(boundary_orth_ECA,2,1); %place boundaries matrix after boundaries matrix for making 'a loop' while choosing the right diameter
        boundary_search = boundary_enlarged(index,:); %subset of boundary matrix, only values relevant for search
        distances = point_to_line_distance(boundary_search, boundaryj, centroid); %distances from line boundaryj-centroid to the points of boundary_search 

        [val,ind] = min(distances); %find value and index of boundarypoint closest to line 
        diameter(j) = pdist([boundaryj;boundary_search(ind,:)]); %Calculate diameter using eucledian distance:
        diameter_plot_ECA{i}(j,:) = [boundaryj boundary_search(ind,:)]; %Save coordinates of diameter
    end
    
    %regionprops:
    stats_ECA = regionprops(binary_ECA,'Area','Circularity','Eccentricity','EquivDiameter','MajorAxisLength','MinorAxisLength','Orientation','Centroid');
    
    diameters_ECA{i} = diameter; %save all diameter values per slice
    mindiameters_ECA(i) = min(diameters_ECA{i}); %smallest diameter for each slice
    maxdiameters_ECA(i) = max(diameters_ECA{i}); %largest diameter for each slice
    avgdiameters_ECA(i) = mean(diameters_ECA{i}); %average diameter for each slice
    area_ECA(i) = pi*((0.5*avgdiameters_ECA(i))^2); %area of vessel for each slice
         
    s_ECA(i)= stats_ECA;
end

%% Pre processing oblique slices

%aantal slices (n_slices) is gelijk laagste slicenummer waar ICA of ECA
%gekozen is. 
if isempty(sliceNr_ICA) && isempty(sliceNr_ECA)
    n_slices = 50; %--------------------------------------------------------DIT NOG INVULLEN!!!
elseif isempty(sliceNr_ICA)
    n_slices = sliceNr_ECA;
elseif isempty(sliceNr_ECA)
    n_slices = sliceNr_ICA;
elseif sliceNr_ICA > sliceNr_ECA
    n_slices = sliceNr_ECA;
elseif sliceNr_ECA > sliceNr_ICA
    n_slices = sliceNr_ICA; 
elseif sliceNr_ICA == sliceNr_ECA
    n_slices = sliceNr_ICA;
end

%Crop center array to smallest amount of slices: 
center_ICA = center_ICA(1:n_slices,:);
center_ECA = center_ECA(1:n_slices,:);

% Make large CTA for oblique slices
size_largeCTA = 10;
lCTA = originalCTA(:,:,lowerbound-size_largeCTA:lowerbound+n_slices+size_largeCTA); %largeCTA
lCTA = flip(lCTA,3); %flip 3rd dimension of CTA such that slice 1 is the lowest slice (containing CCA)

%% Create oblique slices on ORIGINAL dataset for ICA
%Change 3rd column of centerpoints for largeCTA by adding <size_largeCTA> slices:
lcenter_ICA = center_ICA;
lcenter_ICA = [lcenter_ICA(:,1:2),lcenter_ICA(:,3)+size_largeCTA]; %centerpoints are X slices up as well

% Centerpoints/centerline is center_ICA 
% Centerpoints/centerline in largeCTA is lcenter_ICA

% Find slope of centerline ICA: 
for i=1:(size(center_ICA,1)-1)
    d_cent_ICA = lcenter_ICA(i+1,:) - lcenter_ICA(i,:); %centerpoint of next slice minus centerpoint of current slice
    d_center_ICA(i,:) = d_cent_ICA; %store slope in matrix
end

oblique_center_ICA = [lcenter_ICA(:,1) lcenter_ICA(:,2) (lcenter_ICA(:,3)+0.5)]; %take z-coordinate in between slices

% Create oblique slices of ORIGINAL large dataset: 
for i = 1:size(d_center_ICA,1)
    [oblique_slice_ICA,x,y,z] = obliqueslice(lCTA,oblique_center_ICA(i,:),d_center_ICA(i,:), ...
                    'Method','nearest');%extract a slice from CTA
    oblique_slices_ICA{i} = oblique_slice_ICA;
    xen{i} = x;
    yen{i} = y;
    zen{i} = z;
    
    % Find index number in 2D slice of the 3D centerpoint-coordinate:
    pointLinearIndexInSlice = find( ceil(x)==ceil(oblique_center_ICA(i,1)) & ...
                                    ceil(y)==ceil(oblique_center_ICA(i,2)) & ...
                                    ceil(z)==ceil(oblique_center_ICA(i,3)));
                                
   % If there is no solution, look for a coordinate close to the
   % centerpoint:
   if isempty(pointLinearIndexInSlice)
       pointLinearIndexInSlice = find(  (ceil(x)==ceil(oblique_center_ICA(i,1))|ceil(x)==(ceil(oblique_center_ICA(i,1))-1)|ceil(x)==(ceil(oblique_center_ICA(i,1))+1)) &...
                                        (ceil(y)==ceil(oblique_center_ICA(i,2))|ceil(y)==(ceil(oblique_center_ICA(i,2))-1)|ceil(y)==(ceil(oblique_center_ICA(i,2))+1)) &...
                                        (ceil(z)==ceil(oblique_center_ICA(i,3))|ceil(z)==(ceil(oblique_center_ICA(i,3))-1)|ceil(z)==(ceil(oblique_center_ICA(i,3))+1)));
   end
   
   % Define the row and column coordinate in the slice (go from index to
   % coordinate)
   for j = 1:size(pointLinearIndexInSlice,1) %there could be more than one index values
       [pointRow(j),pointColumn(j)] = ind2sub(size(oblique_slice_ICA),pointLinearIndexInSlice(j));
       if size(pointRow,2)>1            %if there are more than 1 index values
           pointRow = pointRow(1);      %choose the first one
           pointColumn = pointColumn(1);
       end
   end
   
   % Save variable in matrices: 
   pointIndex(i) = pointLinearIndexInSlice(1);
   ptRow_ICA(i) = pointRow;
   ptColumn_ICA(i) = pointColumn;           
end


%% Create oblique slices on ORIGINAL dataset for ECA
%Change 3rd column of centerpoints for largeCTA by adding <size_largeCTA> slices:
lcenter_ECA = center_ECA;
lcenter_ECA = [lcenter_ECA(:,1:2),lcenter_ECA(:,3)+size_largeCTA]; %centerpoints are X slices up as well

% Centerpoints/centerline is center_ECA 
% Centerpoints/centerline in largeCTA is lcenter_ICA

% Find slope of centerline ECA: 
for i=1:(size(center_ECA,1)-1)
    d_cent_ECA = lcenter_ECA(i+1,:) - lcenter_ECA(i,:); %centerpoint of next slice minus centerpoint of current slice
    d_center_ECA(i,:) = d_cent_ECA; %store slope in matrix
end

oblique_center_ECA = [lcenter_ECA(:,1) lcenter_ECA(:,2) (lcenter_ECA(:,3)+0.5)]; %take z-coordinate in between slices

% Create oblique slices of ORIGINAL dataset: 
for i = 1:size(d_center_ECA,1)
    [oblique_slice_ECA,x,y,z] = obliqueslice(lCTA,oblique_center_ECA(i,:),d_center_ECA(i,:), ...
                    'Method','nearest');%extract a slice from CTA
    oblique_slices_ECA{i} = oblique_slice_ECA;
    xen{i} = x;
    yen{i} = y;
    zen{i} = z;
    
    % Find index number in 2D slice of the 3D centerpoint-coordinate:
    pointLinearIndexInSlice = find( ceil(x)==ceil(oblique_center_ECA(i,1)) & ...
                                    ceil(y)==ceil(oblique_center_ECA(i,2)) & ...
                                    ceil(z)==ceil(oblique_center_ECA(i,3)));
                                
   % If there is no solution, look for a coordinate close to the
   % centerpoint:
   if isempty(pointLinearIndexInSlice)
       pointLinearIndexInSlice = find(  (ceil(x)==ceil(oblique_center_ECA(i,1))|ceil(x)==(ceil(oblique_center_ECA(i,1))-1)|ceil(x)==(ceil(oblique_center_ECA(i,1))+1)) &...
                                        (ceil(y)==ceil(oblique_center_ECA(i,2))|ceil(y)==(ceil(oblique_center_ECA(i,2))-1)|ceil(y)==(ceil(oblique_center_ECA(i,2))+1)) &...
                                        (ceil(z)==ceil(oblique_center_ECA(i,3))|ceil(z)==(ceil(oblique_center_ECA(i,3))-1)|ceil(z)==(ceil(oblique_center_ECA(i,3))+1)));
   end
   
   % Define the row and column coordinate in the slice (go from index to
   % coordinate)
   for j = 1:size(pointLinearIndexInSlice,1) %there could be more than one index values
       [pointRow(j),pointColumn(j)] = ind2sub(size(oblique_slice_ECA),pointLinearIndexInSlice(j));
       if size(pointRow,2)>1            %if there are more than 1 index values
           pointRow = pointRow(1);      %choose the first one
           pointColumn = pointColumn(1);
       end
   end
   
   % Save variable in matrices: 
   pointIndex(i) = pointLinearIndexInSlice(1);
   ptRow_ECA(i) = pointRow;
   ptColumn_ECA(i) = pointColumn;           
end
    
    
%% Segmentation on oblique slices for ICA
% oblique_slices_ICA = oblique slices van ICA
% ptRow_ICA and ptColumn_ICA zijn centerpoints. Deze verschillen dus van
% oblique_center_ICA in dat ze gelden voor de oblique slice, en niet voor
% het 3D beeld. 
obliqueslice_center_ICA = [ptColumn_ICA' ptRow_ICA']; %centerpoints for oblique slice in 1 array

for i=1:size(oblique_slices_ICA,2)
    current_oblique_slice = oblique_slices_ICA{i};
    binary_oblique_slice = (current_oblique_slice >= lumenMin_app) & ...
                            (current_oblique_slice <= lumenMax_app); % create binary of current slice
    binary_oblique_slices_ICA{i} = binary_oblique_slice; %save binaries in matrix
    %Smooth: blur image and threshold at 0.5
    width = 3;
    kernel = ones(width) / width^2;
    blurryImage_oblique = conv2(double(binary_oblique_slice), kernel, 'same');
    smoothImage_oblique = blurryImage_oblique > 0.5;
    smoothImages_oblique{i} = smoothImage_oblique;
    
    [B,L] = bwboundaries(smoothImage_oblique); %find boundaries and labels of all regions

    % Find the label number where centerpoint is located:
    L_ICA = L(obliqueslice_center_ICA(i,2),obliqueslice_center_ICA(i,1));
    % If label == 0, then find nearest white pixel and give
    % binary_oblique_ICA that label:
    if L_ICA==0
        %calculate distance to all white pixels and choose pixels with
        %smallest distance        
        [x,y] = find(smoothImage_oblique==1); % find all nonzeropixels
        nonzeropixels = [x,y]; % array with non-zero pixels
        distances = sqrt((x - obliqueslice_center_ICA(i,2)).^2 + ...
                            (y - obliqueslice_center_ICA(i,1)).^2); % distances to nonzero pixels:
        [~,sdidx] = min(distances); %linear index of pixel with smallest distance 
        sd = nonzeropixels(sdidx,:); %row-column index of pixel with smallest distance
        L_ICA2 = L(sd(1),sd(2)); % label value of the pixel with smallest distance
        binary_oblique_ICA = L==L_ICA2; %binary with only the label where pixel with smallest distance is located 
    else  
        % Get binary with only the label where centerpoint is located:
        binary_oblique_ICA = L==L_ICA;
    end
    binaries_oblique_ICA{i} = binary_oblique_ICA; %save binary oblique images of ICA in structure
    
    % Get boundaries from ICA on oblique slices:
    B = bwboundaries(binary_oblique_ICA);
    B = B{1};
    boundaries_oblique_ICA{i} = B; 
    
    % Update centerpoint:
    s = regionprops(binary_oblique_ICA,'centroid');
    centerpoint = cat(1,s.Centroid);
    centerpoints_oblique_ICA(i,:) = centerpoint;
    
%     % Segmentation calcium:
%     % Steps: 
%     % 1. Make rectangle ROI around centerpoint (middle of vessel)
%     % 2. Select ROI on oblique slice (non-segmented)
%     % 3. Threshold within the ROI
%     % 4. Calculate boundaries of calcium
%     
%     %Step1:
%     width_ROI = 30; %circle had radius of 15 so width and height of 30 must be sufficient
%     height_ROI = 30;
%     xco = fix(centerpoint(1)-(0.5*width_ROI)); %corner of ROI - use fix to get integers
%     yco = fix(centerpoint(2)-(0.5*height_ROI)); %other corner of ROI - use fix to get integers
%     
%     %Step2:
%     zeros_matrix = zeros(size(current_oblique_slice)); %zeros matrix with size of oblique slice
%         zeros_matrix(yco:yco+height_ROI,xco:xco+width_ROI)=1; %fill zeros matrix with ones in ROI
%     oblique_slice_ROI = current_oblique_slice.*zeros_matrix; %ROI filled with CTA image
%     
%     %Step3:
%     oblique_slice_calc = (oblique_slice_ROI>calciumMin); %threshold
%     
%     %Step4: 
%     B_calc = bwboundaries(oblique_slice_calc); % find boundaries
%     boundaries_oblique_calc_ICA{i} = B_calc; %save boundary per iteration in structure
end
    
    
%% Segmentation on oblique slices for ECA
% oblique_slices_ECA = oblique slices van ECA
% ptRow_ECA and ptColumn_ECA zijn centerpoints. Deze verschillen dus van
% oblique_center_ECA in dat ze gelden voor de oblique slice, en niet voor
% het 3D beeld. 
obliqueslice_center_ECA = [ptColumn_ECA' ptRow_ECA']; %centerpoints for oblique slice in 1 array

for i=1:size(oblique_slices_ECA,2)
    current_oblique_slice = oblique_slices_ECA{i};
    binary_oblique_slice = (current_oblique_slice >= lumenMin_app) & ...
                            (current_oblique_slice <= lumenMax_app); % create binary of current slice
    binary_oblique_slices_ECA{i} = binary_oblique_slice; %save binaries in matrix
    %Smooth: blur image and threshold at 0.5
    width = 3;
    kernel = ones(width) / width^2;
    blurryImage_oblique = conv2(double(binary_oblique_slice), kernel, 'same');
    smoothImage_oblique = blurryImage_oblique > 0.5;
    smoothImages_oblique{i} = smoothImage_oblique;
    
    [B,L] = bwboundaries(smoothImage_oblique); %find boundaries and labels of all regions

    % Find the label number where centerpoint is located:
    L_ECA = L(obliqueslice_center_ECA(i,2),obliqueslice_center_ECA(i,1));
    % If label == 0, then find nearest white pixel and give
    % binary_oblique_ECA that label:
    if L_ECA==0
        %calculate distance to all white pixels and choose pixels with
        %smallest distance        
        [x,y] = find(smoothImage_oblique==1); % find all nonzeropixels
        nonzeropixels = [x,y]; % array with non-zero pixels
        distances = sqrt((x - obliqueslice_center_ECA(i,2)).^2 + ...
                            (y - obliqueslice_center_ECA(i,1)).^2); % distances to nonzero pixels:
        [~,sdidx] = min(distances); %linear index of pixel with smallest distance 
        sd = nonzeropixels(sdidx,:); %row-column index of pixel with smallest distance
        L_ECA2 = L(sd(1),sd(2)); % label value of the pixel with smallest distance
        binary_oblique_ECA = L==L_ECA2; %binary with only the label where pixel with smallest distance is located 
    else  
        % Get binary with only the label where centerpoint is located:
        binary_oblique_ECA = L==L_ECA;
    end
    binaries_oblique_ECA{i} = binary_oblique_ECA; %save binary oblique images of ECA in structure
    
    % Get boundaries from ECA on oblique slices:
    B = bwboundaries(binary_oblique_ECA);
    B = B{1};
    boundaries_oblique_ECA{i} = B; 
    
    % Update centerpoint:
    s = regionprops(binary_oblique_ECA,'centroid');
    centerpoint = cat(1,s.Centroid);
    centerpoints_oblique_ECA(i,:) = centerpoint;
    
%     % Segmentation calcium:
%     % Steps: 
%     % 1. Make rectangle ROI around centerpoint (middle of vessel)
%     % 2. Select ROI on oblique slice (non-segmented)
%     % 3. Threshold within the ROI
%     % 4. Calculate boundaries of calcium
%     
%     %Step1:
%     width_ROI = 30; %circle had radius of 15 so width and height of 30 must be sufficient
%     height_ROI = 30;
%     xco = fix(centerpoint(1)-(0.5*width_ROI)); %corner of ROI - use fix to get integers
%     yco = fix(centerpoint(2)-(0.5*height_ROI)); %other corner of ROI - use fix to get integers
%     
%     %Step2:
%     zeros_matrix = zeros(size(current_oblique_slice)); %zeros matrix with size of oblique slice
%     zeros_matrix(yco:yco+height_ROI,xco:xco+width_ROI)=1; %fill zeros matrix with ones in ROI
%     oblique_slice_ROI = current_oblique_slice.*zeros_matrix; %ROI filled with CTA image
%     
%     %Step3:
%     oblique_slice_calc = (oblique_slice_ROI>calciumMin); %threshold
%     
%     %Step4: 
%     B_calc = bwboundaries(oblique_slice_calc); % find boundaries
%     boundaries_oblique_calc_ECA{i} = B_calc; %save boundary per iteration in structure
end    
    
%% Calculate diameters ICA and other properties
     
for i=1:length(boundaries_oblique_ICA)
    boundary_oblique_ICA = boundaries_oblique_ICA{i}; %select current boundary
    obliqueSlice_ICA = oblique_slices_ICA{i}; %select current oblique slice
    binaryObliqueSlice_ICA = binaries_oblique_ICA{i}; %select current binary oblique slice
    centroid_oblique_ICA = [centerpoints_oblique_ICA(i,2) centerpoints_oblique_ICA(i,1)]; %coordinates of current centroid
    
    % For each boundary point find distance to boundary point at opposite
    % side and choose closest boundary point for calculation of diameter:
    clear diameter    
    for j=1:(0.5*length(boundary_oblique_ICA))
        boundaryj = boundary_oblique_ICA(j,:); %current boundarypoint

        %search for closest point on boundary:
        index = j+floor(0.25*length(boundary_oblique_ICA)):(j+floor(0.25*length(boundary_oblique_ICA))+floor(0.5*length(boundary_oblique_ICA))); %indexes of search area
        boundary_enlarged = repmat(boundary_oblique_ICA,2,1); %place boundaries matrix after boundaries matrix for making 'a loop' while choosing the right diameter
        boundary_search = boundary_enlarged(index,:); %subset of boundary matrix, only values relevant for search
        distances = point_to_line_distance(boundary_search, boundaryj, centroid_oblique_ICA); %distances from line boundaryj-centroid to the points of boundary_search 

        [val,ind] = min(distances); %find value and index of boundarypoint closest to line 
        diameter(j) = pdist([boundaryj;boundary_search(ind,:)]); %Calculate diameter using eucledian distance:
        diameter_plot_oblique_ICA{i}(j,:) = [boundaryj boundary_search(ind,:)]; %Save coordinates of diameter
    end
    
    %regionprops:
    stats_ICA_oblique = regionprops(binaryObliqueSlice_ICA,'Area','Circularity','Eccentricity','EquivDiameter','MajorAxisLength','MinorAxisLength','Orientation','Centroid');
    
    diameters_ICA_oblique{i} = diameter; %save all diameter values per slice
    mindiameters_ICA_oblique(i) = min(diameters_ICA_oblique{i}); %smallest diameter for each slice
    maxdiameters_ICA_oblique(i) = max(diameters_ICA_oblique{i}); %largest diameter for each slice
    avgdiameters_ICA_oblique(i) = mean(diameters_ICA_oblique{i}); %average diameter for each slice
    area_ICA_oblique(i) = pi*((0.5*avgdiameters_ICA_oblique(i))^2); %area of vessel for each slice
    s_ICA_oblique(i) = stats_ICA_oblique;
end
    
    
%% Calculate diameters ECA
     
for i=1:length(boundaries_oblique_ECA)
    boundary_oblique_ECA = boundaries_oblique_ECA{i}; %select current boundary
    obliqueSlice_ECA = oblique_slices_ECA{i}; %select current oblique slice
    binaryObliqueSlice_ECA = binaries_oblique_ECA{i}; %select current binary oblique slice
    centroid_oblique_ECA = [centerpoints_oblique_ECA(i,2) centerpoints_oblique_ECA(i,1)]; %coordinates of current centroid
    
    % For each boundary point find distance to boundary point at opposite
    % side and choose closest boundary point for calculation of diameter:
    clear diameter    
    for j=1:(0.5*length(boundary_oblique_ECA))
        boundaryj = boundary_oblique_ECA(j,:); %current boundarypoint

        %search for closest point on boundary:
        index = j+floor(0.25*length(boundary_oblique_ECA)):(j+floor(0.25*length(boundary_oblique_ECA))+floor(0.5*length(boundary_oblique_ECA))); %indexes of search area
        boundary_enlarged = repmat(boundary_oblique_ECA,2,1); %place boundaries matrix after boundaries matrix for making 'a loop' while choosing the right diameter
        boundary_search = boundary_enlarged(index,:); %subset of boundary matrix, only values relevant for search
        distances = point_to_line_distance(boundary_search, boundaryj, centroid_oblique_ECA); %distances from line boundaryj-centroid to the points of boundary_search 

        [val,ind] = min(distances); %find value and index of boundarypoint closest to line 
        diameter(j) = pdist([boundaryj;boundary_search(ind,:)]); %Calculate diameter using eucledian distance:
        diameter_plot_oblique_ECA{i}(j,:) = [boundaryj boundary_search(ind,:)]; %Save coordinates of diameter
    end
    
    %regionprops:
    stats_ECA_oblique = regionprops(binaryObliqueSlice_ECA,'Area','Circularity','Eccentricity','EquivDiameter','MajorAxisLength','MinorAxisLength','Orientation','Centroid');
    
    diameters_ECA_oblique{i} = diameter; %save all diameter values per slice
    mindiameters_ECA_oblique(i) = min(diameters_ECA_oblique{i}); %smallest diameter for each slice
    maxdiameters_ECA_oblique(i) = max(diameters_ECA_oblique{i}); %largest diameter for each slice
    avgdiameters_ECA_oblique(i) = mean(diameters_ECA_oblique{i}); %average diameter for each slice
    area_ECA_oblique(i) = pi*((0.5*avgdiameters_ECA_oblique(i))^2); %area of vessel for each slice 
    s_ECA_oblique(i) = stats_ECA_oblique;
end


%% Save variables in workspace
% Ask for left or right:
prompt = {'Scannumber?','Left or right carotid artery?'};
answer = inputdlg(prompt);
%Save variables in .mat file to reopen in app_showSmallestArea
filename = sprintf('20240115_Segmentation3D_carotis_test_%s_%s_%d_%d',char(answer(1)),char(answer(2)),lumenMin_app,lumenMax_app);
filename_mat = [filename '.mat'];
savedir = 'C:\Users\Astrid\surfdrive\Documents\MATLAB\Carotis segmentatie\SegmentationResults\';
tic 
save([savedir filename_mat])
toc    
%% Results 
% % Select lower and upper slice of bifurcation area
% % To be used for selecting area for smallest diameter
% app_running = app_selectBifurcation;%launch app
% waitfor(app_running)  % wait for app to close
% % output of app: roi_pos or roi_pos2: a 4-element vector with position of roi [xmin
% % ymin width height]
% 
% if exist('roi_pos2') && ~isempty(roi_pos2) %if roi_pos2 does exist and is not empty
%     roi_pos = roi_pos2; 
% end
% % Bifurcation selection is done in originalCTA. Conversion for CTA:
% lowerslice_bif = fix(roi_pos(2))-lowerbound;
% upperslice_bif = fix(roi_pos(2)+roi_pos(4))-lowerbound;
% 
% % mindiameters_ICA_bif = mindiameters_ICA(lowerslice_bif:upperslice_bif); % mindiameters_ICA_bif are all minimal diameters in bifurcation region
% % [~,mindiameter_ICA_idx] = find(mindiameters_ICA_bif==min(mindiameters_ICA_bif(mindiameters_ICA_bif>1))); % find slice with smallest diameter in bifurcation region that is not zero
% % sliceNr_mindiameter_ICA = mindiameter_ICA_idx+lowerslice_bif-1; %convert slice number to CTA slice number
% 
% results_matrix = [mindiameters_ICA' maxdiameters_ICA' avgdiameters_ICA' area_ICA'];   
% results_table = struct2table(s_ICA);
% results_table.minD = mindiameters_ICA';
% results_table.maxD = maxdiameters_ICA';
% results_table.areaBound = area_ICA'; 
% results_table.sliceNr = (1:length(boundaries_oblique_ICA))';
% 
% t_results = results_table(:,{'sliceNr','minD','maxD','EquivDiameter','MinorAxisLength','MajorAxisLength','areaBound','Area','Circularity','Eccentricity'});
% 
% t_results_bif = t_results(lowerslice_bif:upperslice_bif,:);
% 

%% Save results
% % Ask for left or right:
% prompt = {'Scannumber?','Left or right carotid artery?'};
% answer = inputdlg(prompt);
% %Save variables in .mat file to reopen in app_showSmallestArea
% filename = sprintf('%s_%s',char(answer(1)),char(answer(2)));
% filename_mat = [filename '.mat'];
% filename_xlsx = [filename '.xlsx'];
% savedir = 'C:\Users\Astrid\surfdrive\Documents\MATLAB\Carotis segmentatie\SegmentationResults\';
% save([savedir filename_mat], ...
%     'CTA','originalCTA','lowerbound','upperbound','center_ICA','boundaries_ICA','boundaries_oblique_ICA',...
%     'oblique_slices_ICA','centerpoints_oblique_ICA','diameters_ICA','diameter_plot_ICA','t_results',...
%     't_results_bif','lowerslice_bif','upperslice_bif','lumenMin_app','lumenMax_app','x_CCA','y_CCA',...
%     'x_ICA','y_ICA','x_ECA','y_ECA')
% 
% %Save tables as xlsx file
% writetable(t_results,[savedir 't_results_' filename_xlsx ])
% writetable(t_results_bif,[savedir 't_results_bif_' filename_xlsx ])