clear all; clc; %#ok<CLALL>
DATA = load('data.mat');

DATA.ROK_UR(DATA.ROK_UR == 9997) = NaN;

DATA.WTPfpU(DATA.WTPfpU == 401) = Inf;
DATA.WTPhpU(DATA.WTPhpU == 401) = Inf;

INPUT.bounds = [DATA.WTPfpL DATA.WTPfpU];
% 1-400 €  month or less, 2-600 - 800 €  month, 3-800 - 1000 €  month, 4-1000 - 1600 € month, 5-1600  € month or more

DATA.inc_norm = DATA.INCOME;
DATA.inc_norm(~isnan(DATA.INCOME)) = (DATA.INCOME(~isnan(DATA.INCOME)) - mean(DATA.INCOME(~isnan(DATA.INCOME)))) ./ std(DATA.INCOME(~isnan(DATA.INCOME)));
% INPUT.X = [2015-DATA.ROK_UR, DATA.EDU==2, DATA.EDU==3,DATA.EDU==4,DATA.inc_norm];
% INPUT.NamesX = {'age','edu = 2','edu = 3','edu = 4','income normalized'};
% INPUT.X = [2015-DATA.ROK_UR, DATA.inc_norm];
% INPUT.NamesX = {'age','income normalized'};

INPUT.Spike = 1; % allow for jump density at 0

INPUT.HessEstFix = 2; % available options: 0 - retained from optimization, 1 - BHHH based, 2 - numerical high-precision Jacobian based, 3 - numerical Hessian based
INPUT.SimStats = 1; % simulate descriptive statistics of the distribution

% INPUT.WT = rand(size(INPUT.bounds,1),1); % weight observations

Distributions = {...
    0  'Normal'; ...
    1  'Logistic'; ...
    2  'Extreme_Value'; ...
    3  'Generalized_Extreme_Value'; ...
    4  'tLocationScale'; ...
    5  'Uniform'; ...
    6  'Johnson_SU'; ...
    
    10  'Expotential'; ...
    11  'Lognormal'; ...
    12  'Loglogistic'; ...
    13  'Weibull'; ...
    14  'Rayleigh'; ...
    15  'Gamma'; ...
    16  'BirnbaumSaunders'; ...
    17  'Generalized_Pareto'; ...
    18  'Inverse_Gaussian'; ...
    19  'Nakagami'; ...
    20  'Rician'; ...
    21  'Johnson_SB'; ...
    22  'Johnson_SL'; ...
    
    31  'Poisson'; ...
    32  'Negative_Binomial'
    };


% this is for testing:
INPUT.bounds(2029,2) = Inf;
INPUT.bounds(2028,2) = 0;
INPUT.bounds(2027,1) = -Inf;

for i = 1:size(Distributions,1);
    Results.(Distributions{i,2}) = DistFit(INPUT,Distributions{i,1});
end

for i = 1:size(Distributions,1);
    LL(i,1:2) = {Distributions{i,2}, Results.(Distributions{i,2}).fval}; %#ok<SAGROW>
end

% If the Excel file with the same name exists, its first worksheet will
% be cleared.
if exist(fullfile(pwd,['DistFit',Distributions{i,2},'.xls'])) == 2;
    Excel = actxserver('Excel.Application');
    Workbook = Excel.Workbooks.Open(fullfile(pwd,['DistFit',Distributions{i,2},'.xls']));
    Excel.ActiveWorkBook.Sheets.Item(1).Cells.Clear;
    Workbook.Save;
    Excel.Workbook.Close;
    invoke(Excel, 'Quit');
    delete(Excel)
end

xlswrite(['DistFit',Distributions{i,2},'.xls'], Results.(Distributions{i,2}).R_out)

Excel = actxserver('Excel.Application');
Workbook = Excel.Workbooks.Open(fullfile(pwd,['DistFit',Distributions{i,2},'.xls']));
Worksheet = Workbook.Sheets.Item(1);

Worksheet.Columns.Item(1).columnWidth = 14;
Worksheet.Columns.Item(3).columnWidth = 3;
Worksheet.Columns.Item(7).columnWidth = 3;
Worksheet.Columns.Item(11).columnWidth = 3;

% hExcel.Cells.Select;
% hExcel.Cells.EntireColumn.AutoFit;

Worksheet.Activate;
ActiveWorksheetRange = Excel.Activesheet.get('Range','A:Z');
ActiveWorksheetRange.NumberFormat = '0,0000'; % If Excel uses dots as decimal places, this should be adjusted here.

Workbook.Save
Workbook.Close
Excel.Quit


% col = 1+4*(1*any(i == [10,14,31]) + 2*any(i == [0:2,5,11:13,15,16,18:21,32]) + 3*any(i == [3,4,17,22]) + 4*any(i == [6]))+4*INPUT.Spike;
% if col == 5
%     range = 'A:E';
% elseif col == 9
%     range = 'A:I';
% elseif col == 13
%     range = 'A:M';
% elseif col == 17
%     range = 'A:Q';
% elseif col == 21
%     range = 'A:U';
% end
%
% % xlsAutoFitCol(['DistFit',Distributions{i,2},'.xls'], 'Arkusz1',range)
% xlsAutoFitCol2(['DistFit',Distributions{i,2},'.xls'], range)