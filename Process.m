function [all_results] = Process(path_to_TIL_folder, type, path_to_excel, blurFilterType, s, a, NBack, imgcut, NLvl, open_close)
theFiles = dir(fullfile(append(path_to_TIL_folder,type), '*.png'));
num2 = length(theFiles);

excel_table = readtable(path_to_excel, 'VariableNamingRule', 'preserve');

struct_ness = cell(2,1);
names = cell(2,1);
age = cell(2,1);
race = cell(2,1);
tumor_stage = cell(2,1);
status = cell(2,1);
os = cell(2,1);
os_time = cell(2,1);
DSS = cell(2,1);
DSS_time = cell(2,1);
DFI = cell(2,1);
DFI_time = cell(2,1);
PFI = cell(2,1);
PFI_time = cell(2,1);
patches_num = cell(2,1);
tll_perc = cell(2,1);
counter = 1;
if NLvl == 2
    for o = 1:num2
        baseFileName = theFiles(o).name;
        fullFileName = fullfile(theFiles(o).folder, baseFileName);
        img = imread(fullFileName);
        img2 = double(img(:,:,1));

        if strcmp(open_close, 'open')
            mask = img(:,:,1) == 255 & img(:,:,2)==255 & img(:,:,3)==255;
            img2(mask) = NaN;
            se = strel('square', 2);
            img_negative = 255 - img2;
            img2 = imopen(img_negative, se);
            img2 = 255 - img2;
            img2(mask) = 255;
        elseif strcmp(open_close, 'close')
            mask = img(:,:,1) == 255 & img(:,:,2)==255 & img(:,:,3)==255;
            img2(mask) = NaN;
            se = strel('square', 2);
            img_negative = 255 - img2;
            img2 = imclose(img_negative, se);
            img2 = 255 - img2;
            img2(mask) = 255;
        else
        end

        if strcmp(blurFilterType, 'gaussian')
            mask = img2 == 255;
            img2(mask) = 127.5;
            img2 = imgaussfilt(img2, s);
            img2(mask) = 255;
        elseif strcmp(blurFilterType, 'average')
            mask = img2 == 255;
            img2(mask) = 127.5;
            h = fspecial('average', [a, a]);
            img2 = imfilter(img2, h);
            img2(mask) = 255;
        elseif strcmp(blurFilterType, 'median')
            mask = img2 == 255;
            img2(mask) = 127.5;
            img2 = medfilt2(img2, [2, 2]);
            img2(mask) = 255;
        elseif strcmp(blurFilterType, 'geomean')
            mask = img2 == 255;
            img2(mask) = 127.5;
            img2 = colfilt(img2, [2 2], 'sliding', @geomean);
            img2(mask) = 255;
        end

        if strcmp(NBack,  'y')
            mask = img(:,:,1) == 255 & img(:,:,2)==255 & img(:,:,3)==255;
            img2(mask) = NaN;
        end


        index = find(string(excel_table.bcr_patient_barcode(:)) == extractBetween(baseFileName,1,12));
        tabulate_img2 = tabulate(img2(:));

        if (isempty(index) == 1)
            bigboystr = append(baseFileName,' NOT FOUND skiping...');
            disp(bigboystr);
            continue
        else

            try
                if strcmp(imgcut, 'y')
                    [h, w] = size(img2);
                    h_half = floor(h / 2);
                    w_half = floor(w / 2);

                    im1 = img2(1:h_half, 1:w_half);
                    im2 = img2(1:h_half, w_half+1:end);
                    im3 = img2(h_half+1:end, 1:w_half);
                    im4 = img2(h_half+1:end, w_half+1:end);
                    glcm1 = graycomatrix(im1, 'NumLevels', NLvl, 'Offset', [0 1; -1 1; -1 0; -1 -1], 'Symmetric', false, 'GrayLimits', []);
                    glcm2 = graycomatrix(im2, 'NumLevels', NLvl, 'Offset', [0 1; -1 1; -1 0; -1 -1], 'Symmetric', false, 'GrayLimits', []);
                    glcm3 = graycomatrix(im3, 'NumLevels', NLvl, 'Offset', [0 1; -1 1; -1 0; -1 -1], 'Symmetric', false, 'GrayLimits', []);
                    glcm4 = graycomatrix(im4, 'NumLevels', NLvl, 'Offset', [0 1; -1 1; -1 0; -1 -1], 'Symmetric', false, 'GrayLimits', []);

                    gcm = (glcm1 + glcm2 + glcm3 + glcm4) / 4;
                else
                    gcm = graycomatrix(img2,'NumLevels',NLvl,'Offset',[0 1; -1 1; -1 0; -1 -1],'Symmetric',false,'GrayLimits',[]);
                end
                gcm1 = gcm(:,:,1);
                gcm2 = gcm(:,:,2);
                gcm3 = gcm(:,:,3);
                gcm4 = gcm(:,:,4);

                struct_ness{counter,1} = gcm1(1,1) + gcm2(1,1) + gcm3(1,1) + gcm4(1,1);
                struct_ness{counter,2} = gcm1(2,2) + gcm2(2,2) + gcm3(2,2) + gcm4(2,2);

                names{counter} = baseFileName;
                age{counter} = excel_table.age_at_initial_pathologic_diagnosis(index);
                race{counter} = excel_table.race(index);
                tumor_stage{counter} = excel_table.ajcc_pathologic_tumor_stage(index);
                status{counter} = excel_table.vital_status(index);
                os{counter} = excel_table.OS(index);
                os_time{counter} = excel_table.OS_time(index);
                DSS{counter} = excel_table.DSS(index);
                DSS_time{counter} = excel_table.DSS_time(index);
                DFI{counter} = excel_table.DFI(index);
                DFI_time{counter} = excel_table.DFI_time(index);
                PFI{counter} = excel_table.PFI(index);
                PFI_time{counter} = excel_table.PFI_time(index);
                tll_perc{counter} = tabulate_img2(2,3);
                patches_num{counter} = tabulate_img2(2,2);

                counter = counter + 1;

            catch
                fprintf(2,append(baseFileName,' FAILED TO CALCULATE GLCM skiping...\n'));
                continue
            end
        end
    end


    tab = horzcat(struct_ness);
    tab(:,1) = num2cell(cell2mat(tab(:,1))/max(cell2mat(tab(:,1))));
    tab(:,2) = num2cell(cell2mat(tab(:,2))/max(cell2mat(tab(:,2))));
    tab(:,1) = num2cell(cell2mat(tab(:,1))+cell2mat(tab(:,2)));

    co_occurance_sorted_table = splitvars(sortrows(table(names, tab, age, race, tumor_stage, ...
        status,os, os_time, DSS, DSS_time, DFI, DFI_time, PFI, PFI_time, ...
        patches_num, tll_perc), -2), 'tab', 'NewVariableNames', {'struct_ness_1', 'struct_ness_2'});


    % Przechowanie wyników
    all_results = co_occurance_sorted_table;
elseif NLvl == 4
    for o = 1:num2
        baseFileName = theFiles(o).name;
        fullFileName = fullfile(theFiles(o).folder, baseFileName);
        img = imread(fullFileName);
        img2 = double(img(:,:,1));

        if strcmp(open_close, 'open')
            mask = img(:,:,1) == 255 & img(:,:,2)==255 & img(:,:,3)==255;
            img2(mask) = NaN;
            se = strel('square', 2);
            img_negative = 255 - img2;
            img2 = imopen(img_negative, se);
            img2 = 255 - img2;
            img2(mask) = 255;
        elseif strcmp(open_close, 'close')
            mask = img(:,:,1) == 255 & img(:,:,2)==255 & img(:,:,3)==255;
            img2(mask) = NaN;
            se = strel('square', 2);
            img_negative = 255 - img2;
            img2 = imclose(img_negative, se);
            img2 = 255 - img2;
            img2(mask) = 255;
        else
        end

        if strcmp(blurFilterType, 'gaussian')
            mask = img2 == 255;
            img2(mask) = 127.5;
            img2 = imgaussfilt(img2, s);
            img2(mask) = 255;
        elseif strcmp(blurFilterType, 'average')
            mask = img2 == 255;
            img2(mask) = 127.5;
            h = fspecial('average', [a, a]);
            img2 = imfilter(img2, h);
            img2(mask) = 255;
        elseif strcmp(blurFilterType, 'median')
            mask = img2 == 255;
            img2(mask) = 127.5;
            img2 = medfilt2(img2, [2, 2]);
            img2(mask) = 255;
        elseif strcmp(blurFilterType, 'geomean')
            mask = img2 == 255;
            img2(mask) = 127.5;
            img2 = colfilt(img2, [2 2], 'sliding', @geomean);
            img2(mask) = 255;
        end

        if strcmp(NBack,  'y')
            mask = img(:,:,1) == 255 & img(:,:,2)==255 & img(:,:,3)==255;
            img2(mask) = NaN;
        end


        index = find(string(excel_table.bcr_patient_barcode(:)) == extractBetween(baseFileName,1,12));
        tabulate_img2 = tabulate(img2(:));

        if (isempty(index) == 1)
            bigboystr = append(baseFileName,' NOT FOUND skiping...');
            disp(bigboystr);
            continue
        else

            try
                if strcmp(imgcut, 'y')
                    [h, w] = size(img2);
                    h_half = floor(h / 2);
                    w_half = floor(w / 2);

                    im1 = img2(1:h_half, 1:w_half);
                    im2 = img2(1:h_half, w_half+1:end);
                    im3 = img2(h_half+1:end, 1:w_half);
                    im4 = img2(h_half+1:end, w_half+1:end);
                    glcm1 = graycomatrix(im1, 'NumLevels', NLvl, 'Offset', [0 1; -1 1; -1 0; -1 -1], 'Symmetric', false, 'GrayLimits', []);
                    glcm2 = graycomatrix(im2, 'NumLevels', NLvl, 'Offset', [0 1; -1 1; -1 0; -1 -1], 'Symmetric', false, 'GrayLimits', []);
                    glcm3 = graycomatrix(im3, 'NumLevels', NLvl, 'Offset', [0 1; -1 1; -1 0; -1 -1], 'Symmetric', false, 'GrayLimits', []);
                    glcm4 = graycomatrix(im4, 'NumLevels', NLvl, 'Offset', [0 1; -1 1; -1 0; -1 -1], 'Symmetric', false, 'GrayLimits', []);

                    gcm = (glcm1 + glcm2 + glcm3 + glcm4) / 4;
                else
                    gcm = graycomatrix(img2,'NumLevels',NLvl,'Offset',[0 1; -1 1; -1 0; -1 -1],'Symmetric',false,'GrayLimits',[]);
                end
                M = [2 1 0 0
                    1 1 0 0
                    0 0 1 1
                    0 0 1 2];
                gcm1 = gcm(:,:,1) .* M;
                gcm2 = gcm(:,:,2) .* M;
                gcm3 = gcm(:,:,3) .* M;
                gcm4 = gcm(:,:,4) .* M;

                struct_ness{counter,1} = (gcm1(1,1) + gcm1(2,1) + gcm1(1,2) + gcm1(2,2)) + (gcm2(1,1) + gcm2(2,1) + gcm2(1,2) + gcm2(2,2)) + (gcm3(1,1) + gcm3(2,1) + gcm3(1,2) + gcm3(2,2)) + (gcm4(1,1) + gcm4(2,1) + gcm4(1,2) + gcm4(2,2));
                struct_ness{counter,2} = (gcm1(4,4) + gcm1(3,4) + gcm1(4,3) + gcm1(3,3)) + (gcm2(4,4) + gcm2(3,4) + gcm2(4,3) + gcm2(3,3)) + (gcm3(4,4) + gcm3(3,4) + gcm3(4,3) + gcm3(3,3)) + (gcm4(4,4) + gcm4(3,4) + gcm4(4,3) + gcm4(3,3));

                names{counter} = baseFileName;
                age{counter} = excel_table.age_at_initial_pathologic_diagnosis(index);
                race{counter} = excel_table.race(index);
                tumor_stage{counter} = excel_table.ajcc_pathologic_tumor_stage(index);
                status{counter} = excel_table.vital_status(index);
                os{counter} = excel_table.OS(index);
                os_time{counter} = excel_table.OS_time(index);
                DSS{counter} = excel_table.DSS(index);
                DSS_time{counter} = excel_table.DSS_time(index);
                DFI{counter} = excel_table.DFI(index);
                DFI_time{counter} = excel_table.DFI_time(index);
                PFI{counter} = excel_table.PFI(index);
                PFI_time{counter} = excel_table.PFI_time(index);
                tll_perc{counter} = tabulate_img2(2,3);
                patches_num{counter} = tabulate_img2(2,2);

                counter = counter + 1;

            catch
                fprintf(2,append(baseFileName,' FAILED TO CALCULATE GLCM skiping...\n'));
                continue
            end
        end
    end


    tab = horzcat(struct_ness);
    tab(:,1) = num2cell(cell2mat(tab(:,1))/max(cell2mat(tab(:,1))));
    tab(:,2) = num2cell(cell2mat(tab(:,2))/max(cell2mat(tab(:,2))));
    tab(:,1) = num2cell(cell2mat(tab(:,1))+cell2mat(tab(:,2)));

    co_occurance_sorted_table = splitvars(sortrows(table(names, tab, age, race, tumor_stage, ...
        status,os, os_time, DSS, DSS_time, DFI, DFI_time, PFI, PFI_time, ...
        patches_num, tll_perc), -2), 'tab', 'NewVariableNames', {'struct_ness_1', 'struct_ness_2'});


    % Przechowanie wyników
    all_results = co_occurance_sorted_table;
elseif NLvl == 6
    for o = 1:num2
        baseFileName = theFiles(o).name;
        fullFileName = fullfile(theFiles(o).folder, baseFileName);
        img = imread(fullFileName);
        img2 = double(img(:,:,1));

        if strcmp(open_close, 'open')
            mask = img(:,:,1) == 255 & img(:,:,2)==255 & img(:,:,3)==255;
            img2(mask) = NaN;
            se = strel('square', 2);
            img_negative = 255 - img2;
            img2 = imopen(img_negative, se);
            img2 = 255 - img2;
            img2(mask) = 255;
        elseif strcmp(open_close, 'close')
            mask = img(:,:,1) == 255 & img(:,:,2)==255 & img(:,:,3)==255;
            img2(mask) = NaN;
            se = strel('square', 2);
            img_negative = 255 - img2;
            img2 = imclose(img_negative, se);
            img2 = 255 - img2;
            img2(mask) = 255;
        else
        end

        if strcmp(blurFilterType, 'gaussian')
            mask = img2 == 255;
            img2(mask) = 127.5;
            img2 = imgaussfilt(img2, s);
            img2(mask) = 255;
        elseif strcmp(blurFilterType, 'average')
            mask = img2 == 255;
            img2(mask) = 127.5;
            h = fspecial('average', [a, a]);
            img2 = imfilter(img2, h);
            img2(mask) = 255;
        elseif strcmp(blurFilterType, 'median')
            mask = img2 == 255;
            img2(mask) = 127.5;
            img2 = medfilt2(img2, [2, 2]);
            img2(mask) = 255;
        elseif strcmp(blurFilterType, 'geomean')
            mask = img2 == 255;
            img2(mask) = 127.5;
            img2 = colfilt(img2, [2 2], 'sliding', @geomean);
            img2(mask) = 255;
        end

        if strcmp(NBack,  'y')
            mask = img(:,:,1) == 255 & img(:,:,2)==255 & img(:,:,3)==255;
            img2(mask) = NaN;
        end


        index = find(string(excel_table.bcr_patient_barcode(:)) == extractBetween(baseFileName,1,12));
        tabulate_img2 = tabulate(img2(:));

        if (isempty(index) == 1)
            bigboystr = append(baseFileName,' NOT FOUND skiping...');
            disp(bigboystr);
            continue
        else

            try
                if strcmp(imgcut, 'y')
                    [h, w] = size(img2);
                    h_half = floor(h / 2);
                    w_half = floor(w / 2);

                    im1 = img2(1:h_half, 1:w_half);
                    im2 = img2(1:h_half, w_half+1:end);
                    im3 = img2(h_half+1:end, 1:w_half);
                    im4 = img2(h_half+1:end, w_half+1:end);
                    glcm1 = graycomatrix(im1, 'NumLevels', NLvl, 'Offset', [0 1; -1 1; -1 0; -1 -1], 'Symmetric', false, 'GrayLimits', []);
                    glcm2 = graycomatrix(im2, 'NumLevels', NLvl, 'Offset', [0 1; -1 1; -1 0; -1 -1], 'Symmetric', false, 'GrayLimits', []);
                    glcm3 = graycomatrix(im3, 'NumLevels', NLvl, 'Offset', [0 1; -1 1; -1 0; -1 -1], 'Symmetric', false, 'GrayLimits', []);
                    glcm4 = graycomatrix(im4, 'NumLevels', NLvl, 'Offset', [0 1; -1 1; -1 0; -1 -1], 'Symmetric', false, 'GrayLimits', []);

                    gcm = (glcm1 + glcm2 + glcm3 + glcm4) / 4;
                else
                    gcm = graycomatrix(img2,'NumLevels',NLvl,'Offset',[0 1; -1 1; -1 0; -1 -1],'Symmetric',false,'GrayLimits',[]);
                end
                M = [4 2 1 0 0 0
                    2 2 1 0 0 0
                    1 1 1 0 0 0
                    0 0 0 1 1 1
                    0 0 0 1 2 2
                    0 0 0 1 2 4];
                gcm1 = gcm(:,:,1) .* M;
                gcm2 = gcm(:,:,2) .* M;
                gcm3 = gcm(:,:,3) .* M;
                gcm4 = gcm(:,:,4) .* M;

                struct_ness{counter,1} = (gcm1(1,1) + gcm1(2,1) + gcm1(3,1) + gcm1(1,2) + gcm1(2,2) + gcm1(3,2) + gcm1(1,3) + gcm1(2,3) + gcm1(3,3)) + (gcm2(1,1) + gcm2(2,1) + gcm2(3,1) + gcm2(1,2) + gcm2(2,2) + gcm2(3,2) + gcm2(1,3) + gcm2(2,3) + gcm2(3,3)) + (gcm3(1,1) + gcm3(2,1) + gcm3(3,1) + gcm3(1,2) + gcm3(2,2) + gcm3(3,2) + gcm3(1,3) + gcm3(2,3) + gcm3(3,3)) + (gcm4(1,1) + gcm4(2,1) + gcm4(3,1) + gcm4(1,2) + gcm4(2,2) + gcm4(3,2) + gcm4(1,3) + gcm4(2,3) + gcm4(3,3));
                struct_ness{counter,2} = (gcm1(6,6) + gcm1(5,6) + gcm1(4,6) + gcm1(6,5) + gcm1(5,5) + gcm1(4,5) + gcm1(6,4) + gcm1(5,4) + gcm1(4,4)) + (gcm2(6,6) + gcm2(5,6) + gcm2(4,6) + gcm2(6,5) + gcm2(5,5) + gcm2(4,5) + gcm2(6,4) + gcm2(5,4) + gcm2(4,4)) + (gcm3(6,6) + gcm3(5,6) + gcm3(4,6) + gcm3(6,5) + gcm3(5,5) + gcm3(4,5) + gcm3(6,4) + gcm3(5,4) + gcm3(4,4)) + (gcm4(6,6) + gcm4(5,6) + gcm4(4,6) + gcm4(6,5) + gcm4(5,5) + gcm4(4,5) + gcm4(6,4) + gcm4(5,4) + gcm4(4,4));

                names{counter} = baseFileName;
                age{counter} = excel_table.age_at_initial_pathologic_diagnosis(index);
                race{counter} = excel_table.race(index);
                tumor_stage{counter} = excel_table.ajcc_pathologic_tumor_stage(index);
                status{counter} = excel_table.vital_status(index);
                os{counter} = excel_table.OS(index);
                os_time{counter} = excel_table.OS_time(index);
                DSS{counter} = excel_table.DSS(index);
                DSS_time{counter} = excel_table.DSS_time(index);
                DFI{counter} = excel_table.DFI(index);
                DFI_time{counter} = excel_table.DFI_time(index);
                PFI{counter} = excel_table.PFI(index);
                PFI_time{counter} = excel_table.PFI_time(index);
                tll_perc{counter} = tabulate_img2(2,3);
                patches_num{counter} = tabulate_img2(2,2);

                counter = counter + 1;

            catch
                fprintf(2,append(baseFileName,' FAILED TO CALCULATE GLCM skiping...\n'));
                continue
            end
        end
    end

    tab = horzcat(struct_ness);
    tab(:,1) = num2cell(cell2mat(tab(:,1))/max(cell2mat(tab(:,1))));
    tab(:,2) = num2cell(cell2mat(tab(:,2))/max(cell2mat(tab(:,2))));
    tab(:,1) = num2cell(cell2mat(tab(:,1))+cell2mat(tab(:,2)));

    co_occurance_sorted_table = splitvars(sortrows(table(names, tab, age, race, tumor_stage, ...
        status,os, os_time, DSS, DSS_time, DFI, DFI_time, PFI, PFI_time, ...
        patches_num, tll_perc), -2), 'tab', 'NewVariableNames', {'struct_ness_1', 'struct_ness_2'});


    % Przechowanie wyników
    all_results = co_occurance_sorted_table;
elseif NLvl == 8
    for o = 1:num2
        baseFileName = theFiles(o).name;
        fullFileName = fullfile(theFiles(o).folder, baseFileName);
        img = imread(fullFileName);
        img2 = double(img(:,:,1));

        if strcmp(open_close, 'open')
            mask = img(:,:,1) == 255 & img(:,:,2)==255 & img(:,:,3)==255;
            img2(mask) = NaN;
            se = strel('square', 2);
            img_negative = 255 - img2;
            img2 = imopen(img_negative, se);
            img2 = 255 - img2;
            img2(mask) = 255;
        elseif strcmp(open_close, 'close')
            mask = img(:,:,1) == 255 & img(:,:,2)==255 & img(:,:,3)==255;
            img2(mask) = NaN;
            se = strel('square', 2);
            img_negative = 255 - img2;
            img2 = imclose(img_negative, se);
            img2 = 255 - img2;
            img2(mask) = 255;
        else
        end

        if strcmp(blurFilterType, 'gaussian')
            mask = img2 == 255;
            img2(mask) = 127.5;
            img2 = imgaussfilt(img2, s);
            img2(mask) = 255;
        elseif strcmp(blurFilterType, 'average')
            mask = img2 == 255;
            img2(mask) = 127.5;
            h = fspecial('average', [a, a]);
            img2 = imfilter(img2, h);
            img2(mask) = 255;
        elseif strcmp(blurFilterType, 'median')
            mask = img2 == 255;
            img2(mask) = 127.5;
            img2 = medfilt2(img2, [2, 2]);
            img2(mask) = 255;
        elseif strcmp(blurFilterType, 'geomean')
            mask = img2 == 255;
            img2(mask) = 127.5;
            img2 = colfilt(img2, [2 2], 'sliding', @geomean);
            img2(mask) = 255;
        end

        if strcmp(NBack,  'y')
            mask = img(:,:,1) == 255 & img(:,:,2)==255 & img(:,:,3)==255;
            img2(mask) = NaN;
        end


        index = find(string(excel_table.bcr_patient_barcode(:)) == extractBetween(baseFileName,1,12));
        tabulate_img2 = tabulate(img2(:));

        if (isempty(index) == 1)
            bigboystr = append(baseFileName,' NOT FOUND skiping...');
            disp(bigboystr);
            continue
        else

            try
                if strcmp(imgcut, 'y')
                    [h, w] = size(img2);
                    h_half = floor(h / 2);
                    w_half = floor(w / 2);

                    im1 = img2(1:h_half, 1:w_half);
                    im2 = img2(1:h_half, w_half+1:end);
                    im3 = img2(h_half+1:end, 1:w_half);
                    im4 = img2(h_half+1:end, w_half+1:end);
                    glcm1 = graycomatrix(im1, 'NumLevels', NLvl, 'Offset', [0 1; -1 1; -1 0; -1 -1], 'Symmetric', false, 'GrayLimits', []);
                    glcm2 = graycomatrix(im2, 'NumLevels', NLvl, 'Offset', [0 1; -1 1; -1 0; -1 -1], 'Symmetric', false, 'GrayLimits', []);
                    glcm3 = graycomatrix(im3, 'NumLevels', NLvl, 'Offset', [0 1; -1 1; -1 0; -1 -1], 'Symmetric', false, 'GrayLimits', []);
                    glcm4 = graycomatrix(im4, 'NumLevels', NLvl, 'Offset', [0 1; -1 1; -1 0; -1 -1], 'Symmetric', false, 'GrayLimits', []);

                    gcm = (glcm1 + glcm2 + glcm3 + glcm4) / 4;
                else
                    gcm = graycomatrix(img2,'NumLevels',NLvl,'Offset',[0 1; -1 1; -1 0; -1 -1],'Symmetric',false,'GrayLimits',[]);
                end
                M = [4 2 1 0 0 0 0 0
                    2 2 1 0 0 0 0 0
                    1 1 1 0 0 0 0 0
                    0 0 0 0 0 0 0 0
                    0 0 0 0 0 0 0 0
                    0 0 0 0 0 1 1 1
                    0 0 0 0 0 1 2 2
                    0 0 0 0 0 1 2 4];
                gcm1 = gcm(:,:,1) .* M;
                gcm2 = gcm(:,:,2) .* M;
                gcm3 = gcm(:,:,3) .* M;
                gcm4 = gcm(:,:,4) .* M;

                struct_ness{counter,1} = (gcm1(1,1) + gcm1(2,1) + gcm1(3,1) + gcm1(1,2) + gcm1(2,2) + gcm1(3,2) + gcm1(1,3) + gcm1(2,3) + gcm1(3,3)) + (gcm2(1,1) + gcm2(2,1) + gcm2(3,1) + gcm2(1,2) + gcm2(2,2) + gcm2(3,2) + gcm2(1,3) + gcm2(2,3) + gcm2(3,3)) + (gcm3(1,1) + gcm3(2,1) + gcm3(3,1) + gcm3(1,2) + gcm3(2,2) + gcm3(3,2) + gcm3(1,3) + gcm3(2,3) + gcm3(3,3)) + (gcm4(1,1) + gcm4(2,1) + gcm4(3,1) + gcm4(1,2) + gcm4(2,2) + gcm4(3,2) + gcm4(1,3) + gcm4(2,3) + gcm4(3,3));
                struct_ness{counter,2} = (gcm1(8,8) + gcm1(7,8) + gcm1(6,8) + gcm1(8,7) + gcm1(7,7) + gcm1(6,7) + gcm1(8,6) + gcm1(7,6) + gcm1(6,6)) + (gcm2(8,8) + gcm2(7,8) + gcm2(6,8) + gcm2(8,7) + gcm2(7,7) + gcm2(6,7) + gcm2(8,6) + gcm2(7,6) + gcm2(6,6)) + (gcm3(8,8) + gcm3(7,8) + gcm3(6,8) + gcm3(8,7) + gcm3(7,7) + gcm3(6,7) + gcm3(8,6) + gcm3(7,6) + gcm3(6,6)) + (gcm4(8,8) + gcm4(7,8) + gcm4(6,8) + gcm4(8,7) + gcm4(7,7) + gcm4(6,7) + gcm4(8,6) + gcm4(7,6) + gcm4(6,6));


                names{counter} = baseFileName;
                age{counter} = excel_table.age_at_initial_pathologic_diagnosis(index);
                race{counter} = excel_table.race(index);
                tumor_stage{counter} = excel_table.ajcc_pathologic_tumor_stage(index);
                status{counter} = excel_table.vital_status(index);
                os{counter} = excel_table.OS(index);
                os_time{counter} = excel_table.OS_time(index);
                DSS{counter} = excel_table.DSS(index);
                DSS_time{counter} = excel_table.DSS_time(index);
                DFI{counter} = excel_table.DFI(index);
                DFI_time{counter} = excel_table.DFI_time(index);
                PFI{counter} = excel_table.PFI(index);
                PFI_time{counter} = excel_table.PFI_time(index);
                tll_perc{counter} = tabulate_img2(2,3);
                patches_num{counter} = tabulate_img2(2,2);

                counter = counter + 1;

            catch
                fprintf(2,append(baseFileName,' FAILED TO CALCULATE GLCM skiping...\n'));
                continue
            end
        end
    end

    tab = horzcat(struct_ness);
    tab(:,1) = num2cell(cell2mat(tab(:,1))/max(cell2mat(tab(:,1))));
    tab(:,2) = num2cell(cell2mat(tab(:,2))/max(cell2mat(tab(:,2))));
    tab(:,1) = num2cell(cell2mat(tab(:,1))+cell2mat(tab(:,2)));

    co_occurance_sorted_table = splitvars(sortrows(table(names, tab, age, race, tumor_stage, ...
        status,os, os_time, DSS, DSS_time, DFI, DFI_time, PFI, PFI_time, ...
        patches_num, tll_perc), -2), 'tab', 'NewVariableNames', {'struct_ness_1', 'struct_ness_2'});


    % Przechowanie wyników
    all_results = co_occurance_sorted_table;
end
end