%% open file
fid = 'TCGACRC_expression.xlsx'; % File Name of Excel here
par = detectImportOptions(fid);
par.VariableNames = {'feature'};  % Modify , Narrow search area
genes = readtable(fid,par); genes=table2cell(genes); % 'cell' Output
%% Data Prep
genearray = string(genes);i=0; clc,check =1;
while check; i=i+1;
genename =input('Gene to Compare: ','s'); rownum = find(genearray== genename); % = 3414
if isempty(rownum);disp('Gene cannot be found or you said you are done');check = 0; 
else disp([genename,' is at ',num2str(rownum),' row. ']);
targenes{i}=genename; end 
% Import gene expression row data
;fid=fopen('TCGACRC_expression.tsv');frewind(fid); %fopen file
for count = 1:rownum, fgetl(fid);end,arr2=fgetl(fid);
subcall=@(A,st) A(st:end);gene_exp =subcall(strsplit(arr2),2); %Data ouput cell
geneexp{i} = str2double(gene_exp); end,loop = i;  %Data output double
geneexp(i:end) = []; targenes(i:end)=[]; % remove entry
for a = 1:numel(geneexp);expvalues(a,:)=geneexp{a};end
%% Clustergram
clustergram(expvalues(:,1:end)')%,'RowLabels',targenes,'ColumnLabels',1:length(expvalues))
% 1. Input expression 2. genes must be cell array of string
input('Press enter to terminate program'), close all;