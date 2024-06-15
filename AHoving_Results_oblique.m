
%% ICA - Create table for saving:
results_matrix_ICA_oblique = [mindiameters_ICA_oblique' maxdiameters_ICA_oblique' avgdiameters_ICA_oblique' area_ICA_oblique'];   
results_table_ICA_oblique = struct2table(s_ICA_oblique);
results_table_ICA_oblique.minD = mindiameters_ICA_oblique';
results_table_ICA_oblique.maxD = maxdiameters_ICA_oblique';
results_table_ICA_oblique.avgD = avgdiameters_ICA_oblique';
results_table_ICA_oblique.areaBound = area_ICA_oblique'; 
results_table_ICA_oblique.sliceNr = (1:length(boundaries_oblique_ICA))';

t_results_ICA_oblique = results_table_ICA_oblique(:,{'sliceNr','minD','maxD','avgD','EquivDiameter','MinorAxisLength','MajorAxisLength','areaBound','Area','Circularity','Eccentricity'});

%% ECA - Create table for saving:
results_matrix_ECA_oblique = [mindiameters_ECA_oblique' maxdiameters_ECA_oblique' avgdiameters_ECA_oblique' area_ECA_oblique'];   
results_table_ECA_oblique = struct2table(s_ECA_oblique);
results_table_ECA_oblique.minD = mindiameters_ECA_oblique';
results_table_ECA_oblique.maxD = maxdiameters_ECA_oblique';
results_table_ECA_oblique.avgD = avgdiameters_ECA_oblique';
results_table_ECA_oblique.areaBound = area_ECA_oblique'; 
results_table_ECA_oblique.sliceNr = (1:length(boundaries_oblique_ECA))';

t_results_ECA_oblique = results_table_ECA_oblique(:,{'sliceNr','minD','maxD','avgD','EquivDiameter','MinorAxisLength','MajorAxisLength','areaBound','Area','Circularity','Eccentricity'});

%% ICA - Create table for saving mm:
%Convert variables to mm: 
mindiameters_ICA_oblique_mm = mindiameters_ICA_oblique*dicom_information.PixelSpacing(1);
maxdiameters_ICA_oblique_mm = maxdiameters_ICA_oblique*dicom_information.PixelSpacing(1);
avgdiameters_ICA_oblique_mm = avgdiameters_ICA_oblique*dicom_information.PixelSpacing(1);

%Calculate Area_mm from regionprops Area:
Area_ICA_oblique_mm = results_table_ICA_oblique.Area*dicom_information.PixelSpacing(1)*dicom_information.PixelSpacing(2);
Area_ICA_oblique_mm = table(Area_ICA_oblique_mm,'VariableNames',"area");
%results_table_mm = [results_table(:,1) Area_mm results_table(:,2:end)];

%Create table in mm for saving:
results_matrix_ICA_oblique_mm = [mindiameters_ICA_oblique_mm' maxdiameters_ICA_oblique_mm' avgdiameters_ICA_oblique_mm'];   
results_table_ICA_oblique_mm = table;
results_table_ICA_oblique_mm.minD = mindiameters_ICA_oblique_mm';
results_table_ICA_oblique_mm.maxD = maxdiameters_ICA_oblique_mm';
results_table_ICA_oblique_mm.avgD = avgdiameters_ICA_oblique_mm';
%results_table_mm.Area = Area_mm;
results_table_ICA_oblique_mm.sliceNr = (1:length(boundaries_oblique_ICA))';
t_results_ICA_oblique_mm = results_table_ICA_oblique_mm(:,{'sliceNr','minD','maxD','avgD'});
t_results_ICA_oblique_mm = [t_results_ICA_oblique_mm Area_ICA_oblique_mm];


%% ECA - Create table for saving mm
%Convert variables to mm: 
mindiameters_ECA_oblique_mm = mindiameters_ECA_oblique*dicom_information.PixelSpacing(1);
maxdiameters_ECA_oblique_mm = maxdiameters_ECA_oblique*dicom_information.PixelSpacing(1);
avgdiameters_ECA_oblique_mm = avgdiameters_ECA_oblique*dicom_information.PixelSpacing(1);

%Calculate Area_mm from regionprops Area:
Area_ECA_oblique_mm = results_table_ECA_oblique.Area*dicom_information.PixelSpacing(1)*dicom_information.PixelSpacing(2);
Area_ECA_oblique_mm = table(Area_ECA_oblique_mm,'VariableNames',"area");
%results_table_mm = [results_table(:,1) Area_mm results_table(:,2:end)];

%Create table in mm for saving:
results_matrix_ECA_oblique_mm = [mindiameters_ECA_oblique_mm' maxdiameters_ECA_oblique_mm' avgdiameters_ECA_oblique_mm'];   
results_table_ECA_oblique_mm = table;
results_table_ECA_oblique_mm.minD = mindiameters_ECA_oblique_mm';
results_table_ECA_oblique_mm.maxD = maxdiameters_ECA_oblique_mm';
results_table_ECA_oblique_mm.avgD = avgdiameters_ECA_oblique_mm';
%results_table_mm.Area = Area_mm;
results_table_ECA_oblique_mm.sliceNr = (1:length(boundaries_oblique_ECA))';
t_results_ECA_oblique_mm = results_table_ECA_oblique_mm(:,{'sliceNr','minD','maxD','avgD'});
t_results_ECA_oblique_mm = [t_results_ECA_oblique_mm Area_ECA_oblique_mm];



%% Plot slice number vs. diameter and area 
%Pre-processing:
sliceNrs_ICA_oblique = 1:length(boundaries_oblique_ICA); %number of slices ICA
sliceNrs_ECA_oblique = 1:length(boundaries_oblique_ECA); %number of slices ECA

Area_ECA_oblique_mm = table2array(Area_ECA_oblique_mm);
Area_ICA_oblique_mm = table2array(Area_ICA_oblique_mm);

%Find slice number of bifurcation (where areas differ in size)
if length(Area_ECA_oblique_mm)<length(Area_ICA_oblique_mm)
    Area_ECA_oblique_mm_enl = [Area_ECA_oblique_mm' zeros(1,numel(Area_ICA_oblique_mm)-numel(Area_ECA_oblique_mm))]'; % make arrays same size
    diffICA_ECA_oblique = Area_ICA_oblique_mm - Area_ECA_oblique_mm_enl; %difference between area ICA and area ECA
    bifurcation_oblique = find(diffICA_ECA_oblique,1); %index of first value not zero, so this is the slice of bifurcation
elseif length(Area_ICA_oblique_mm)<length(Area_ECA_oblique_mm)
    Area_ICA_oblique_mm_enl = [Area_ICA_oblique_mm' zeros(1,numel(Area_ECA_oblique_mm)-numel(Area_ICA_oblique_mm))]'; % make arrays same size
    diffICA_ECA_oblique = Area_ECA_oblique_mm - Area_ICA_oblique_mm_enl; %difference between area ICA and area ECA
    bifurcation_oblique = find(diffICA_ECA_oblique,1); %index of first value not zero, so this is the slice of bifurcation
else
    diffICA_ECA_oblique = Area_ECA_oblique_mm - Area_ICA_oblique_mm; %difference between area ICA and area ECA
    bifurcation_oblique = find(diffICA_ECA_oblique,1); %index of first value not zero, so this is the slice of bifurcation
end

%Plot avg, min and max diameters and Areas
figure;plot(sliceNrs_ICA_oblique,avgdiameters_ICA_oblique_mm, sliceNrs_ICA_oblique,mindiameters_ICA_oblique_mm, sliceNrs_ICA_oblique,maxdiameters_ICA_oblique_mm, sliceNrs_ICA_oblique,Area_ICA_oblique_mm)
hold on;plot(sliceNrs_ECA_oblique,avgdiameters_ECA_oblique_mm, sliceNrs_ECA_oblique,mindiameters_ECA_oblique_mm, sliceNrs_ECA_oblique,maxdiameters_ECA_oblique_mm, sliceNrs_ECA_oblique,Area_ECA_oblique_mm)
xline(bifurcation_oblique,'-',{'Bifurcation'},'LineWidth',1.5)
xlabel('Slice number')
ylabel('Diameter (mm) / Area (mm2)')
legend('Average ICA','MinD ICA','MaxD ICA','Area ICA','Average ECA','MinD ECA','MaxD ECA','Area ECA')


%% calculate minimal diameter
% if minD=0 AND Area=dicom_information.PixelSpacing(1)*dicom_information.PixelSpacing(2) OR  
% Area=2*dicom_information.PixelSpacing(1)*dicom_information.PixelSpacing(2)=> minD = dicom_information.PixelSpacing(1)

%sten gives slicenumber(s) of stenosis
if sten == 0
    stenrange = [bifurcation_oblique sliceNrs_ICA_oblique(end)];
else
    stenrange = [sten'-5 sten'+5];
end

lm_oblique = islocalmin(Area_ICA_oblique_mm);
for i=1:size(stenrange,1)
    idx_lm_oblique_temp{i} = find(lm_oblique(stenrange(i,1):stenrange(i,2))); %find slice numbers of local minima in stenosis range
    idx_lm_oblique{i} = length(Area_ICA_oblique_mm(1:stenrange(i,1))) + idx_lm_oblique_temp{i} -1;%def slice numbers of local minima
    [min_Area_ICA_oblique(i),idx_min_Area_ICA_oblique_temp] = min(Area_ICA_oblique_mm(idx_lm_oblique{i}));
    idx_min_Area_ICA_oblique(i) = idx_lm_oblique{i}(idx_min_Area_ICA_oblique_temp); %slice nr(s) of smallest Area
    minD_lm_oblique(i) = mindiameters_ICA_oblique_mm(idx_min_Area_ICA_oblique(i)); %min diameter at slice(s) with smallest area
end

if length(minD_lm_oblique)>1
    for i=1:length(minD_lm_oblique)
        area_minD_oblique(i) = Area_ICA_oblique_mm(idx_min_Area_ICA_oblique(i)) + minD_lm_oblique(i); %add area and minDiameter
    end
    [min_area_minD_oblique,idx_minD_oblique_temp] = min(area_minD_oblique); 
    minD_oblique = minD_lm_oblique(idx_minD_oblique_temp); %min diameter at lowest Area+minD
    idx_minD_oblique = idx_min_Area_ICA_oblique(idx_minD_oblique_temp);%slice number min diameter
else
    minD_oblique = minD_lm_oblique; %min diameter
    idx_minD_oblique = idx_min_Area_ICA_oblique; %slicenummer min diameter
end

%if minD==0 give it the value of 1 pixel (in mm):
pixelvalue = dicom_information.PixelSpacing(1)*dicom_information.PixelSpacing(2);
minD_zero=0;
if (minD_oblique == 0) && (Area_ICA_oblique_mm(idx_minD_oblique) == pixelvalue || Area_ICA_oblique_mm(idx_minD_oblique) == 2*pixelvalue)
    minD_oblique = dicom_information.PixelSpacing(1);
    minD_oblique_zero = 1; %to check if this loop is true
end

%% calculate "normal" diameter
%choose largest maxD after lm slice
[maxD_oblique_ref,idx_maxD_oblique]= max(maxdiameters_ICA_oblique_mm(idx_minD_oblique+1:end)); %maxD_ref = referentie diameter
idx_maxD_oblique_ref = idx_maxD_oblique + idx_minD_oblique; %slice number of "normal" diameter

%% calculate stenosegraad
stenosegraad_oblique = 1-minD_oblique/maxD_oblique_ref;

%% Save variables in workspace
% Ask for left or right:
prompt = {'Scannumber?','Left or right carotid artery?'};
answer = inputdlg(prompt);
%Save variables in .mat file to reopen in app_showSmallestArea
filename = sprintf('20240115_Results_oblique_%s_%s_%d_%d',char(answer(1)),char(answer(2)),lumenMin_app,lumenMax_app);
filename_mat = [filename '.mat'];
savedir = 'C:\Users\Astrid\surfdrive\Documents\MATLAB\Carotis segmentatie\SegmentationResults\';
save([savedir filename_mat], 'stenosegraad_oblique','maxD_oblique_ref','idx_maxD_oblique',...
    'minD_oblique','idx_minD_oblique',...
    'sliceNrs_ICA_oblique','avgdiameters_ICA_oblique_mm', 'mindiameters_ICA_oblique_mm', ...
    'maxdiameters_ICA_oblique_mm', 'Area_ICA_oblique_mm',...
    'sliceNrs_ECA_oblique', 'avgdiameters_ECA_oblique_mm', 'mindiameters_ECA_oblique_mm', ...
    'maxdiameters_ECA_oblique_mm', 'Area_ECA_oblique_mm', ...
    't_results_ECA_oblique_mm', 't_results_ICA_oblique_mm', 't_results_ECA_oblique', 't_results_ICA_oblique')






