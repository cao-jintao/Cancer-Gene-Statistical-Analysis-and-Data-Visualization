%% open file
fid = 'TCGACRC_expression.xlsx'; % File Name of Excel here
par = detectImportOptions(fid);
par.VariableNames = {'feature'};  % Modify , Narrow search area
gene = readtable(fid,par);  %Data Out = Table('abc')
%% Import Data % 1. Import Data 2. Name Data 'gene'  3. Import Displayed rows
name='CDH17'; % NAME of target gene % Import name as 'cdh17_exp' and 'gene_exp'
genearray = string(table2cell(gene));genename =input('Gene to Compare: ','s');
numrow = find(genearray == genename); namerow = find(genearray == name); % = 3414
if isempty(numrow);disp('Gene cannot be found');return,end
disp([genename,' is at ',num2str(numrow),' row. ', name,' is at '...
,num2str(namerow),' row.']);%uiimport('TCGACRC_expression.tsv')% File Name(.txt)
% Import gene expression row data
fid=fopen('TCGACRC_expression.tsv');frewind(fid);bool=numrow>namerow;%% fopen FILE
for i = 1:namerow,fgetl(fid);end,arr1=fgetl(fid); frewind(fid);
for i = 1:numrow, fgetl(fid);end,arr2=fgetl(fid);
subcall=@(A,st) A(st:end);cdh17_exp=subcall(strsplit(arr1),2);
gene_exp =subcall(strsplit(arr2),2);
% gene analysis
geneexp = str2double(gene_exp); cdh17exp = str2double(cdh17_exp); 
%%
nsets = 3
for i = 1:nsets-1;eval(sprintf('geneexp_g%d = geneexp(g%d);',i,i));
eval(sprintf('cdh17exp_g%d = cdh17exp(g%d);',i,i));end 
% gene masks (standard deviation)
per95 = @(num) mean(num)-3*std(num);
lval = per95(cdh17exp)<cdh17exp & per95(geneexp)<geneexp;
geneexp=geneexp(lval);cdh17exp=cdh17exp(lval);% low value mask
%%
%{
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
legend(name,'fit',genename,'fit'),title('normalized beta
distribution');hold off
%}
% Correlation between 2 genes
% x_gene1 = hgene(1).YData; no_cdh17 = h17(1).YData; %no. in expression lv.
%%
corr = corrcoef(cdh17exp,geneexp);corr = corr(1,2);
%normcorr = corrcoef(y_cdh17_dist,y_gene_dist);normcorr = normcorr(1,2);
%Stats Calculation (t-test, pvalue)
stats = regstats(cdh17exp,geneexp);
t_stat = stats.tstat.t(2); p_stat = stats.tstat.pval(2);
% Scatter plot + plot fit
figure;scatter(cdh17exp,geneexp,10,linspace(1,10,length(geneexp)),'filled');
axis([0 20 0 20]), ylabel(genename),xlabel(name),hold on
title(strcat(genename,' vs ', name));
deg = 1; %deg = input('Fitting degree(1 for linear): ') % Input
p = polyfit(cdh17exp, geneexp,deg);b = polyval(p,cdh17exp);%POLYNOMIAL FIT
%X = [ones(numel(cdh17exp),1) cdh17exp'];b1= X\y_gene';%linear regression
Rsq2 = 1- sum((geneexp-b).^2)/sum((geneexp-mean(geneexp)).^2); %Coef of R^2
plot(cdh17exp',b,'LineWidth',2),legend('Scatter','fit'),hold off,
textbox2 = sprintf(['Correlation by patient = %.3f\nR^{2} Coef. = %.3f\n'...
    'T Stat: %.2e , p-value: %.2e'],corr,Rsq2,t_stat,p_stat);
annotation('textbox',[.14 .7 .34 .22],'string',textbox2);
%%  Plot graphs of various tumor grades
% Require to load TCGAGRC_clinical.xlsx first for clinicalpathological
% labels
figure;
for i=1:nsets-1; % loop of grades
    subplot(2,3,i); 
    eval(sprintf('cdh17exp=cdh17exp_g%d;geneexp = geneexp_g%d;',i,i))
    scatter(cdh17exp,geneexp,10,linspace(1,10,length(geneexp)),'filled');
    ylabel(genename),xlabel(name),hold on
    title(strcat('grade: ',num2str(sets(i)),',n=',num2str(numel(geneexp)))); 
    deg=1;  stats = regstats(cdh17exp,geneexp);
    corrmat = corrcoef(cdh17exp,geneexp);corr(i) = corrmat(1,2);
    t_stat = stats.tstat.t(2); p_stat = stats.tstat.pval(2);
    p = polyfit(cdh17exp, geneexp,deg);b = polyval(p,cdh17exp);%POLYNOMIAL FIT
    Rsq2 = 1- sum((geneexp-b).^2)/sum((geneexp-mean(geneexp)).^2); %Coef of R^2
    plot(cdh17exp',b,'LineWidth',2),hold off, % ADD LEGEND if neccessary
    text{i} = sprintf(['\nRun number:%d\nCorr = %.3f\nR^{2} coef: %.3f\n' ...
        'TStat: %.3e , p-value: %.3e'],i,corr(i),Rsq2,t_stat,p_stat);
end
boxcoefs = strjoin(string(text),'\n');
textbox2 = sprintf('Grading Classification: %s\n',colname);
figure;annotation('textbox',[0 0 1 1],'string',strcat(textbox2,boxcoefs));
%
for i=1:nsets-1;nrun{i}=sprintf('Grading: %d',sets(i));end;
figure;plot(1:nsets-1,corr,'r-*'); ylabel('correlation');xlabel('runs');
title(strcat('Correlations across Tumor: ',colname,' Grading'));
set(gca,'XTick',1:nsets-1,'XTickLabel',string(nrun));
savegene=input('Press enter to terminate or save gene: ','s');
if isempty(savegene); close all,clc; else,file = fopen('savegene.txt','a+'); 
    fprintf(file,'\n %s %f \n',savegene, corr);fclose(file);
    close all,clc,end
%% z-score scatter group plot
[pc, zscores, pcvars] = pca(expvalues);
figure; scatter(zscores(:,1),zscores(:,2));
xlabel('First Principal Component');ylabel('Second Principal Component');
title('Principal Component Scatter Plot');