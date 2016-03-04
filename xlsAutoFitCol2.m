function xlsAutoFitCol2(filename,varargin)

% set the column to auto fit of the sheet1

%xlsAutoFitCol(filename,sheetname,range)
% Example:
%xlsAutoFitCol('filename','A:F')

options = varargin;

range = varargin{1};
    
[fpath,file,ext] = fileparts(char(filename));

if isempty(fpath)
    fpath = pwd;
end

Excel = actxserver('Excel.Application');
set(Excel,'Visible',0);
Workbook = invoke(Excel.Workbooks, 'open', [fpath filesep file ext]);
invoke(Workbook.Sheets.Item(1),'Activate');

ExAct = Excel.Activesheet;
 
ExActRange = get(ExAct,'Range',range);
ExActRange.Select;

invoke(Excel.Selection.Columns,'Autofit');

invoke(Workbook, 'Save');
invoke(Excel, 'Quit');
delete(Excel);

