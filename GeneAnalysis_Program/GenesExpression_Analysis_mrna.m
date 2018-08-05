%% open file
fid = 'mrna_expression.txt'; %'TCGACRC_expression.xlsx'; % File Name of Excel here
par = detectImportOptions(fid); 
par.VariableNames = {'COMMON'} ;%{'feature'};  % Modify , Narrow search area
gene = readtable(fid,par);  %Data Out = Table('abc')
gene = gene(:,2);
%% Import Data % 1. Import Data 2. Name Data 'gene'  3. Import Displayed rows
name='CDH17'; % NAME of target gene % Import name as 'cdh17_exp' and 'gene_exp'
genearray = string(table2cell(gene));genename =input('Gene to Compare: ','s');
numrow = find(genearray == genename); namerow = find(genearray == name); % = 3414
if length(numrow)>1,numrow = numrow(1);end
if isempty(numrow);disp('Gene cannot be found');return,end
disp([genename,' is at ',num2str(numrow),' row. ', name,' is at '...
,num2str(namerow),' row.']);%uiimport('TCGACRC_expression.tsv')% File Name(.txt)

% Import gene expression row data
fid=fopen('mrna_expression.txt');frewind(fid);bool=numrow>namerow;%% fopen FILE
if bool,a=namerow;b=numrow;else a=numrow;b=namerow;end ;dif=abs(a-b); %create bool
for i = 1:a;fgetl(fid);end;arr1=fgetl(fid);for i=1:dif-1;fgetl(fid);end
arr2=fgetl(fid); subindex=@(A,start) A(start:end);  %## Anonymous function HERE
arr1 = subindex(strsplit(arr1),2);arr2 = subindex(strsplit(arr2),2);
if bool,cdh17_exp=arr1; gene_exp =arr2;else gene_exp=arr1; cdh17_exp=arr2;end

% GENE ANALYSIS
geneexp = str2double(gene_exp); cdh17exp = str2double(cdh17_exp);
% Gene masks
stdscreen = @(num) mean(num)-2*std(num);
nanscreen = isnan(geneexp) & isnan(cdh17exp);
geneexp(nanscreen)=[]; cdh17exp(nanscreen)=[]; %remove nan
lval = stdscreen(geneexp)<geneexp & stdscreen(cdh17exp)<geneexp;
geneexp=geneexp(lval);cdh17exp=cdh17exp(lval);% low value masklow value mask
%
figure; hgene = histfit(geneexp,numel(geneexp),'kernel');     %Figure 1
y_gene_dist = get(hgene(2),'YData');
prop=get(gca,'Children');set(prop(2),'FaceColor',[0 0.5 0.5]); hold on
title(strcat(name,'/',genename,' Expression in CRC patients'))
ylabel('No. of Samples');xlabel('Normalized Expression Level'),legend(genename)
%
h17 = histfit(cdh17exp,numel(cdh17exp),'kernel');%(..,..,'beta','poisson'...etc)
y_cdh17_dist =get(h17(2),'YData'); % distribution y data
legend(genename,'fit',name,'fit'),hold off ;
cdh17peak = max(y_cdh17_dist);genepeak = max(y_gene_dist);
cdh17mean = h17(2).XData(y_cdh17_dist==cdh17peak);
genemean = hgene(2).XData(y_gene_dist==genepeak);
textbox1 = sprintf(['The shift in the expression of %s is '...
    '%.2f %% of %s.\n '],name,cdh17mean/genemean*100,genename);
annotation('textbox',[.14 .7 .34 .2],'string',textbox1)   %text box

% Beta FIGURE
figure;histfit(cdh17exp/max(cdh17exp),numel(cdh17exp),'beta'); 
ylabel('no. of samples');xlabel('expression (normalised)');hold on,
histfit(geneexp/(max(geneexp)),numel(geneexp),'beta');
prop2 = get(gca,'Children'); set(prop2(2),'FaceColor',[0 0.5 0.5]);
legend(name,'fit',genename,'fit'),title('normalized beta distribution');hold off

% Correlation between 2 genes
% x_gene1 = hgene(1).YData; no_cdh17 = h17(1).YData; %no. in expression lv

corr = corrcoef(cdh17exp,geneexp);corr = corr(1,2);
normcorr = corrcoef(y_cdh17_dist,y_gene_dist);normcorr = normcorr(1,2);

% Scatter plot + plot fit
figure;scatter(cdh17exp,geneexp,10,linspace(1,10,length(geneexp)),'filled');
ylabel(genename),xlabel('CDH17'),hold on
title(strcat('Expression plot : ', genename,' vs CDH17'));
deg = 1; %deg = input('Fitting degree(1 for linear): ') % Input
p = polyfit(cdh17exp, geneexp,deg);b = polyval(p,cdh17exp);%POLYNOMIAL FIT
%X = [ones(numel(cdh17exp),1) cdh17exp'];b1= X\y_gene';%linear regression
Rsq2 = 1- sum((geneexp-b).^2)/sum((geneexp-mean(geneexp)).^2); %Coef of R^2
plot(cdh17exp',b,'LineWidth',2),legend('Scatter','fit'),hold off,
textbox2 = sprintf(['Data Correlation by patient = %.3f\n'...
    'The coefficient of determination R^{2} is %.3f'],corr,Rsq2);
annotation('textbox',[.14 .7 .34 .22],'string',textbox2)
input('Press enter to terminate');close all,clc;