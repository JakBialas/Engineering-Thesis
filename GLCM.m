function [data] = GLCM(path_to_TIL_folder, path_to_excel, type)
% Kombinacje Modyfikacji
NLvl_values = [2, 4, 6, 8];%2, 4, 6, 8
NBack_values = {'n', 'y'};%'n','y'
blurFilterType_values = {'n', 'average', 'gaussian', 'median', 'geomean'};%'n', 'average', 'gaussian', 'median', 'geomean'
imgcut_values = {'n', 'y'};%,'n', 'y'
open_close_values = {'n', 'close', 'open'};%'n', 'close', 'open'
s_values = [1, 2, 3];%1, 2, 3
a_values = [3, 5, 7];%3, 5, 7
disp('Beginning GLCM procedure, this may take a while...');

% Inicjacja komórki do przechowania wyniku
all_results = cell(length(NLvl_values) * length(NBack_values) * ...
    length(blurFilterType_values) * length(imgcut_values) * length(s_values) * length(a_values), 1);
result_counter = 1;
% Inicjalizacja komórki do przechowywania opisów
descriptions = cell(length(NLvl_values) * length(NBack_values) * ...
    length(blurFilterType_values) * length(imgcut_values) * length(s_values) * length(a_values), 1);
T = table('Size', [0 1], 'VariableTypes', {'double'}, 'VariableNames', {'Harrells_c_index'});
NT = table('Size', [0 1], 'VariableTypes', {'double'}, 'VariableNames', {'Negatywny_Harrells_c_index'});
P = table('Size', [0 1], 'VariableTypes', {'double'}, 'VariableNames', {'P_values'});
H = table('Size', [0 1], 'VariableTypes', {'double'}, 'VariableNames', {'Hazard_ratio'});
HL = table('Size', [0 1], 'VariableTypes', {'double'}, 'VariableNames', {'Hazard_ratio_in_logarithm'});
C = table('Size', [0 1], 'VariableTypes', {'double'}, 'VariableNames', {'Confidence_intervals'});
CL = table('Size', [0 1], 'VariableTypes', {'double'}, 'VariableNames', {'Confidence_intervals_in_logarithm'});
K = table('Size', [0 1], 'VariableTypes', {'double'}, 'VariableNames', {'Description'});

for i = 1:length(NLvl_values)
    for j = 1:length(NBack_values)
        for k = 1:length(blurFilterType_values)
            for l = 1:length(imgcut_values)
                for p = 1:length(open_close_values)
                    if strcmp(blurFilterType_values(k), 'average')
                        for m = 1:length(a_values)
                            NLvl = NLvl_values(i);
                            NBack = NBack_values{j};
                            blurFilterType = blurFilterType_values{k};
                            imgcut = imgcut_values{l};
                            open_close = open_close_values{p};
                            a = a_values(m);
                            s = 0;
                            disp(['type: ' type ', NLvl: ' num2str(NLvl) ', NBack: ' NBack ', open_close: ' open_close ', blurFilterType: ' blurFilterType ', imgcut: ' imgcut ', a: ' num2str(a)]);
                            % Proces nakładania modyfikacji oraz GLCM
                            all_results{result_counter} = Process(path_to_TIL_folder, type, path_to_excel, blurFilterType, s, a, NBack, imgcut, NLvl, open_close);

                            % Tworzenie opisu
                            description = ['NLvl: ' num2str(NLvl) ', NBack: ' NBack ', open_close: ' open_close ', blurFilterType: ' blurFilterType ', imgcut: ' imgcut ', a: ' num2str(a)];
                            description = cellstr(description);
                            descriptions{result_counter} = description;
                            result_counter = result_counter + 1;
                            disp('Done');
                        end
                    elseif strcmp(blurFilterType_values(k), 'gaussian')
                        for n = 1:length(s_values)
                            NLvl = NLvl_values(i);
                            NBack = NBack_values{j};
                            blurFilterType = blurFilterType_values{k};
                            imgcut = imgcut_values{l};
                            open_close = open_close_values{p};
                            a = 0;
                            s = s_values(n);
                            disp(['type: ' type ', NLvl: ' num2str(NLvl) ', NBack: ' NBack ', open_close: ' open_close ', blurFilterType: ' blurFilterType ', imgcut: ' imgcut ', s: ' num2str(s)]);
                            % Proces nakładania modyfikacji oraz GLCM
                            all_results{result_counter} = Process(path_to_TIL_folder, type, path_to_excel, blurFilterType, s, a, NBack, imgcut, NLvl, open_close);

                            % Tworzenie opisu
                            description = ['NLvl: ' num2str(NLvl) ', NBack: ' NBack ', open_close: ' open_close ', blurFilterType: ' blurFilterType ', imgcut: ' imgcut ', s: ' num2str(s)];
                            description = cellstr(description);
                            descriptions{result_counter} = description;
                            result_counter = result_counter + 1;
                            disp('Done');
                        end
                    elseif strcmp(blurFilterType_values(k), 'median')
                        NLvl = NLvl_values(i);
                        NBack = NBack_values{j};
                        blurFilterType = blurFilterType_values{k};
                        imgcut = imgcut_values{l};
                        open_close = open_close_values{p};
                        a = 0;
                        s = 0;
                        disp(['type: ' type ', NLvl: ' num2str(NLvl) ', NBack: ' NBack ', open_close: ' open_close ', blurFilterType: ' blurFilterType ', imgcut: ' imgcut]);
                        % Proces nakładania modyfikacji oraz GLCM
                        all_results{result_counter} = Process(path_to_TIL_folder, type, path_to_excel, blurFilterType, s, a, NBack, imgcut, NLvl, open_close);

                        % Tworzenie opisu
                        description = ['NLvl: ' num2str(NLvl) ', NBack: ' NBack ', open_close: ' open_close ', blurFilterType: ' blurFilterType ', imgcut: ' imgcut];
                        description = cellstr(description);
                        descriptions{result_counter} = description;
                        result_counter = result_counter + 1;
                        disp('Done');
                    elseif strcmp(blurFilterType_values(k), 'geomean')
                        NLvl = NLvl_values(i);
                        NBack = NBack_values{j};
                        blurFilterType = blurFilterType_values{k};
                        imgcut = imgcut_values{l};
                        open_close = open_close_values{p};
                        a = 0;
                        s = 0;
                        disp(['type: ' type ', NLvl: ' num2str(NLvl) ', NBack: ' NBack ', open_close: ' open_close ', blurFilterType: ' blurFilterType ', imgcut: ' imgcut]);
                        % Proces nakładania modyfikacji oraz GLCM
                        all_results{result_counter} = Process(path_to_TIL_folder, type, path_to_excel, blurFilterType, s, a, NBack, imgcut, NLvl, open_close);

                        % Tworzenie opisu
                        description = ['NLvl: ' num2str(NLvl) ', NBack: ' NBack ', open_close: ' open_close ', blurFilterType: ' blurFilterType ', imgcut: ' imgcut];
                        description = cellstr(description);
                        descriptions{result_counter} = description;
                        result_counter = result_counter + 1;
                        disp('Done');
                    else
                        NLvl = NLvl_values(i);
                        NBack = NBack_values{j};
                        blurFilterType = blurFilterType_values{k};
                        imgcut = imgcut_values{l};
                        open_close = open_close_values{p};
                        a = 0;
                        s = 0;
                        disp(['type: ' type ', NLvl: ' num2str(NLvl) ', NBack: ' NBack ', open_close: ' open_close ', blurFilterType: ' blurFilterType ', imgcut: ' imgcut]);

                        all_results{result_counter} = Process(path_to_TIL_folder, type, path_to_excel, blurFilterType, s, a, NBack, imgcut, NLvl, open_close);

                        % Tworzenie opisu
                        description = ['NLvl: ' num2str(NLvl) ', NBack: ' NBack ', open_close: ' open_close ', blurFilterType: ' blurFilterType ', imgcut: ' imgcut ];
                        description = cellstr(description);
                        descriptions{result_counter} = description;
                        result_counter = result_counter + 1;
                        disp('Done');
                    end
                end
            end
        end
    end
end


% Tworzenie tabeli z wynikami
for i = 1:length(all_results)
    if ~isempty(all_results{i})
        results = all_results{i};
        %Obliczanie wartości analizy przeżycia
        [concordance, minusconcordance, p_values, HR, log_HR, Conf_interv, Conf_interv_log] = Analiza(results, 'struct_ness_1');
        newRow = table(concordance, 'VariableNames', {'Harrells_c_index'});
        newRowN = table(minusconcordance, 'VariableNames', {'Negatywny_Harrells_c_index'});
        newRowp = table(p_values(1,:), 'VariableNames', {'P_values'});
        newRowHR = table(HR, 'VariableNames', {'Hazard_ratio'});
        newRowlog_HR = table(log_HR, 'VariableNames', {'Hazard_ratio_in_logarithm'});
        newRowConf_interv = table(cellstr(num2str(Conf_interv)), 'VariableNames', {'Confidence_intervals'});
        newRowConf_interv_log = table(cellstr(num2str(Conf_interv_log)), 'VariableNames', {'Confidence_intervals_in_logarithm'});
        T = [T; newRow];
        NT = [NT; newRowN];
        P = [P; newRowp];
        H = [H; newRowHR];
        HL = [HL; newRowlog_HR];
        C = [C; newRowConf_interv];
        CL = [CL; newRowConf_interv_log];
    end
end
%Zapisywanie opisu
for i = 1:length(descriptions)
    if ~isempty(descriptions{i})
        results = descriptions{i};
        newRo = table(results, 'VariableNames', {'Description'});
        K = [K; newRo];
    end
end
data = [K, T ,NT, P, H, HL, C, CL];

disp('All combinations processed.');