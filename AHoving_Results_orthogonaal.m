
%% ICA - Create table for saving:
results_matrix_ICA = [mindiameters_ICA' maxdiameters_ICA' avgdiameters_ICA' area_ICA'];   
results_table_ICA = struct2table(s_ICA);
results_table_ICA.minD = mindiameters_ICA';
results_table_ICA.maxD = maxdiameters_ICA';
results_table_ICA.avgD = avgdiameters_ICA';
results_table_ICA.areaBound = area_ICA'; 
results_table_ICA.sliceNr = (1:length(boundaries_ICA))';

t_results_ICA = results_table_ICA(:,{'sliceNr','minD','maxD','avgD','EquivDiameter','MinorAxisLength','MajorAxisLength','areaBound','Area','Circularity','Eccentricity'});

%% ECA - Create table for saving:
results_matrix_ECA = [mindiameters_ECA' maxdiameters_ECA' avgdiameters_ECA' area_ECA'];   
results_table_ECA = struct2table(s_ECA);
results_table_ECA.minD = mindiameters_ECA';
results_table_ECA.maxD = maxdiameters_ECA';
results_table_ECA.avgD = avgdiameters_ECA';
results_table_ECA.areaBound = area_ECA'; 
results_table_ECA.sliceNr = (1:length(boundaries_ECA))';

t_results_ECA = results_table_ECA(:,{'sliceNr','minD','maxD','avgD','EquivDiameter','MinorAxisLength','MajorAxisLength','areaBound','Area','Circularity','Eccentricity'});

%% ICA - Create table for saving mm:
%Convert variables to mm: 
mindiameters_ICA_mm = mindiameters_ICA*dicom_information.PixelSpacing(1);
maxdiameters_ICA_mm = maxdiameters_ICA*dicom_information.PixelSpacing(1);
avgdiameters_ICA_mm = avgdiameters_ICA*dicom_information.PixelSpacing(1);

%Calculate Area_mm from regionprops Area:
Area_ICA_mm = results_table_ICA.Area*dicom_information.PixelSpacing(1)*dicom_information.PixelSpacing(2);
Area_ICA_mm = table(Area_ICA_mm,'VariableNames',"area");
%results_table_mm = [results_table(:,1) Area_mm results_table(:,2:end)];

%Create table in mm for saving:
results_matrix_ICA_mm = [mindiameters_ICA_mm' maxdiameters_ICA_mm' avgdiameters_ICA_mm'];   
results_table_ICA_mm = table;
results_table_ICA_mm.minD = mindiameters_ICA_mm';
results_table_ICA_mm.maxD = maxdiameters_ICA_mm';
results_table_ICA_mm.avgD = avgdiameters_ICA_mm';
%results_table_mm.Area = Area_mm;
results_table_ICA_mm.sliceNr = (1:length(boundaries_ICA))';
t_results_ICA_mm = results_table_ICA_mm(:,{'sliceNr','minD','maxD','avgD'});
t_results_ICA_mm = [t_results_ICA_mm Area_ICA_mm];


%% ECA - Create table for saving mm
%Convert variables to mm: 
mindiameters_ECA_mm = mindiameters_ECA*dicom_information.PixelSpacing(1);
maxdiameters_ECA_mm = maxdiameters_ECA*dicom_information.PixelSpacing(1);
avgdiameters_ECA_mm = avgdiameters_ECA*dicom_information.PixelSpacing(1);

%Calculate Area_mm from regionprops Area:
Area_ECA_mm = results_table_ECA.Area*dicom_information.PixelSpacing(1)*dicom_information.PixelSpacing(2);
Area_ECA_mm = table(Area_ECA_mm,'VariableNames',"area");
%results_table_mm = [results_table(:,1) Area_mm results_table(:,2:end)];

%Create table in mm for saving:
results_matrix_ECA_mm = [mindiameters_ECA_mm' maxdiameters_ECA_mm' avgdiameters_ECA_mm'];   
results_table_ECA_mm = table;
results_table_ECA_mm.minD = mindiameters_ECA_mm';
results_table_ECA_mm.maxD = maxdiameters_ECA_mm';
results_table_ECA_mm.avgD = avgdiameters_ECA_mm';
%results_table_mm.Area = Area_mm;
results_table_ECA_mm.sliceNr = (1:length(boundaries_ECA))';
t_results_ECA_mm = results_table_ECA_mm(:,{'sliceNr','minD','maxD','avgD'});
t_results_ECA_mm = [t_results_ECA_mm Area_ECA_mm];



%% Plot slice number vs. diameter and area 
%Pre-processing:
sliceNrs_ICA = 1:length(boundaries_ICA); %number of slices ICA
sliceNrs_ECA = 1:length(boundaries_ECA); %number of slices ECA

Area_ECA_mm = table2array(Area_ECA_mm);
Area_ICA_mm = table2array(Area_ICA_mm);

%Find slice number of bifurcation (where areas differ in size)
if length(Area_ECA_mm)<length(Area_ICA_mm)
    Area_ECA_mm_enl = [Area_ECA_mm' zeros(1,numel(Area_ICA_mm)-numel(Area_ECA_mm))]'; % make arrays same size
    diffICA_ECA = Area_ICA_mm - Area_ECA_mm_enl; %difference between area ICA and area ECA
    bifurcation = find(diffICA_ECA,1); %index of first value not zero, so this is the slice of bifurcation
elseif length(Area_ICA_mm)<length(Area_ECA_mm)
    Area_ICA_mm_enl = [Area_ICA_mm' zeros(1,numel(Area_ECA_mm)-numel(Area_ICA_mm))]'; % make arrays same size
    diffICA_ECA = Area_ECA_mm - Area_ICA_mm_enl; %difference between area ICA and area ECA
    bifurcation = find(diffICA_ECA,1); %index of first value not zero, so this is the slice of bifurcation
else
    diffICA_ECA = Area_ECA_mm - Area_ICA_mm; %difference between area ICA and area ECA
    bifurcation = find(diffICA_ECA,1); %index of first value not zero, so this is the slice of bifurcation
end

%Plot avg, min and max diameters and Areas
figure;plot(sliceNrs_ICA,avgdiameters_ICA_mm, sliceNrs_ICA,mindiameters_ICA_mm, sliceNrs_ICA,maxdiameters_ICA_mm, sliceNrs_ICA,Area_ICA_mm)
hold on;plot(sliceNrs_ECA,avgdiameters_ECA_mm, sliceNrs_ECA,mindiameters_ECA_mm, sliceNrs_ECA,maxdiameters_ECA_mm, sliceNrs_ECA,Area_ECA_mm)
xline(bifurcation,'-',{'Bifurcation'},'LineWidth',1.5)
xlabel('Slice number')
ylabel('Diameter (mm) / Area (mm2)')
legend('Average ICA','MinD ICA','MaxD ICA','Area ICA','Average ECA','MinD ECA','MaxD ECA','Area ECA')


%% calculate minimal diameter
% if minD=0 AND Area=dicom_information.PixelSpacing(1)*dicom_information.PixelSpacing(2) OR  
% Area=2*dicom_information.PixelSpacing(1)*dicom_information.PixelSpacing(2)=> minD = dicom_information.PixelSpacing(1)

if sten == 0
    stenrange = [bifurcation sliceNrs_ICA(end)];
else
    stenrange = [sten'-5 sten'+5];
end

lm = islocalmin(Area_ICA_mm);
for i=1:size(stenrange,1)
    idx_lm_temp{i} = find(lm(stenrange(i,1):stenrange(i,2))); %find slice numbers of local minima in stenosis range
    idx_lm{i} = length(Area_ICA_mm(1:stenrange(i,1))) + idx_lm_temp{i} -1; %def slice numbers of local minima
    [min_Area_ICA(i),idx_min_Area_ICA_temp] = min(Area_ICA_mm(idx_lm{i}));
    idx_min_Area_ICA(i) = idx_lm{i}(idx_min_Area_ICA_temp); %slice nr(s) of smallest Area
    minD_lm(i) = mindiameters_ICA_mm(idx_min_Area_ICA(i)); %min diameter at slice(s) with smallest area
end

if length(minD_lm)>1
    for i=1:length(minD_lm)
        area_minD(i) = Area_ICA_mm(idx_min_Area_ICA(i)) + minD_lm(i); %add area and minDiameter
    end
    [min_area_minD,idx_minD_temp] = min(area_minD); 
    minD = minD_lm(idx_minD_temp); %min diameter at lowest Area+minD
    idx_minD = idx_min_Area_ICA(idx_minD_temp);%slice number min diameter
else
    minD = minD_lm; %min diameter
    idx_minD = idx_min_Area_ICA; %slicenummer min diameter
end
    
%if minD==0 give it the value of 1 pixel (in mm):
pixelvalue = dicom_information.PixelSpacing(1)*dicom_information.PixelSpacing(2);
minD_zero=0;
if (minD == 0) && (Area_ICA_mm(idx_minD) == pixelvalue || Area_ICA_mm(idx_minD) == 2*pixelvalue)
    minD = dicom_information.PixelSpacing(1);
    minD_zero = 1; %to check if this loop is true
end

%% calculate "normal" diameter
%choose largest maxD after minD slice:
[maxD_ref,idx_maxD]= max(maxdiameters_ICA_mm(idx_minD+1:end)); %maxD_ref = referentie diameter
idx_maxD_ref = idx_maxD + idx_minD; %slice number of "normal" diameter

%% calculate stenosegraad
stenosegraad = 1-minD/maxD_ref;

%% Save variables in workspace
% Ask for left or right:
prompt = {'Scannumber?','Left or right carotid artery?'};
answer = inputdlg(prompt);
%Save variables in .mat file to reopen in app_showSmallestArea
filename = sprintf('20240115_Results_orthogonaal_%s_%s_%d_%d',char(answer(1)),char(answer(2)),lumenMin_app,lumenMax_app);
filename_mat = [filename '.mat'];
savedir = 'C:\Users\Astrid\surfdrive\Documents\MATLAB\Carotis segmentatie\SegmentationResults\';
save([savedir filename_mat], 'stenosegraad','maxD_ref','idx_maxD','minD','idx_minD',...
    'sliceNrs_ICA','avgdiameters_ICA_mm', 'mindiameters_ICA_mm', ...
    'maxdiameters_ICA_mm', 'Area_ICA_mm',...
    'sliceNrs_ECA', 'avgdiameters_ECA_mm', 'mindiameters_ECA_mm', ...
    'maxdiameters_ECA_mm', 'Area_ECA_mm', ...
    't_results_ECA_mm', 't_results_ICA_mm', 't_results_ECA', 't_results_ICA')
    
    






