clear all; close all; warning('off', 'all');
% Ścierzka do folderu z obrazami TIL, należy zamieścić '\' na końcu ścieżki
path_to_TIL_folder = 'C:\Users\jakub\Desktop\DaneTestowe\Dane do testu kodu\';
% Ścierzka do pliku excela zawierająca dane pacjentów
path_to_excel = 'TCGA-CDR-SupplementalTableS1.xlsx';

% Miejsce zapisu pliku excela z wynikami, należy na końcu ścierzki zamieścić '\'
save_results_to = 'C:\Users\jakub\Desktop\DaneTestowe\Dane do testu kodu\';
tic;
% Typy raka wybrane do analizy
typelist = {'Zdjęcia testowe'};
for i = 1:length(typelist)
    type = typelist{i};
    data = GLCM(path_to_TIL_folder, path_to_excel, type);
    data = table(data, 'VariableNames', {'Data'});

    form = [save_results_to, 'Harrells_c_index_data_', type, '.xlsx'];
    data = splitvars(data);
    writetable(data, form);

    disp(['Processed ', type, ' and saved to ', form]);
end

elapsedTimeInSeconds = toc;
% Czas z sekund na minuty
elapsedTimeInMinutes = elapsedTimeInSeconds / 60;

% Czas w minutach
disp(['Czas wykonania kodu: ', num2str(elapsedTimeInMinutes), ' minut']);

disp('All combinations and types processed.');