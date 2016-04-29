function xls = xlsResults(xlsname, content)

% xlsname - the name of the xls document in the format name.xls
% content - what will be saved in the xls document

% If an Excel file with the same name exists, 
% its first worksheet will be cleared.
if exist(fullfile(pwd,xlsname)) == 2;
    Excel = actxserver('Excel.Application');
    Workbook = Excel.Workbooks.Open(fullfile(pwd,xlsname));
    Excel.ActiveWorkBook.Sheets.Item(1).Cells.Clear;
    Workbook.Save;
    Excel.Workbook.Close;
    invoke(Excel, 'Quit');
    delete(Excel)
end

xlswrite(xlsname, content)

Excel = actxserver('Excel.Application');
Workbook = Excel.Workbooks.Open(fullfile(pwd,xlsname));
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