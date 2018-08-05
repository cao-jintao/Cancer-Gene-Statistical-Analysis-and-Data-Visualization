%% (yeast example)
load yeastdata.mat 
fprintf('Amount:%.0f,Test:%s',numel(genes),genes{15})
%% (gene expression set)
% uiimport('TCGACRC_expression.tsv') %manual select
fid = 'TCGACRC_expression.xlsx'; % File Name of Excel here
par = detectImportOptions(fid); 
par.VariableNames = {'feature'};  % Modify , Narrow search area
gene = readtable(fid,par); genearray = string(table2cell(gene)); %
%% (gene expression set)
fid = fopen('TCGACRC_expression.tsv'); i = 1;
geneinput{1} = input('Gene: ','s');
while isempty(geneinput{i})~= 1
i=i+1; geneinput{i}=input('Gene: ','s'); end; geneinput{end} = [];
%% Data prep (geneset)
for m=1:numel(geneinput)-1      %extrating all desired rows
    numrow = find(genearray==geneinput{m})
    for n = 1:numrow;row{m}= fgetl(fid);row{m}=cellstr(strsplit(row{m}));
    end; frewind(fid); 
end
dataset=row{1};
for m = 1:length(row)-1
    dataset = vertcat(dataset,row{m+1});end
dataset = str2double(dataset);
%% Expression vs time (yeast example)
plot(times, expvalues(15,:)),
xlabel('Time (Hours)');ylabel('Log2 Relative Expression Level');
%% ######################## Extract data (gene expression set)
dataset = table2cell(data); % table input code
for i=length(dataset):numel(dataset),if iscategorical(dataset{i}),dataset{i}=[];end,end 
expvalues = cell2mat(dataset(1:end,3:end));genes = dataset(1:end,1);
times = 1:1:length(expvalues);
%% Filtering empty data space
emptySpots = strcmp('EMPTY',genes); %for STRING compare
expvalues(emptySpots,:)=[];genes(emptySpots)=[];
fprintf('\n1st Number: %f',numel(genes))
%
nanIndices = any(isnan(expvalues),2);
expvalues(nanIndices,:)=[]; genes(nanIndices)=[];
fprintf('\n2nd number: %f',numel(genes))
%% High-order MASK
mask = genevarfilter(expvalues);
expvalues=expvalues(mask,:);genes=genes(mask); 
fprintf('\nGene-masked number: %f',numel(genes))
% Low value Filter
[mask, expvalues, genes] = ...
genelowvalfilter(expvalues,genes,'absval',log2(4));
fprintf('\nLowvalue-maksed number: %f',numel(genes))
% Entropy Filter
[mask, expvalues, genes] = ...
geneentropyfilter(expvalues,genes,'prctile',15);
fprintf('\nEntropy-masked number: %f\n',numel(genes))
%% Gene Cluster 
corrDist = pdist(expvalues, 'corr');  % Expression Input Here
clusterTree = linkage(corrDist, 'average'); 
n = 16; clusters = cluster(clusterTree, 'maxclust',n); %Input no. of plots
figure; nr = sqrt(n)% n must be sqare
for c = 1:n
    subplot(nr,nr,c);plot(times,expvalues((clusters == c),:)')
    axis tight
end
suptitle('Hierarchical Clustering of Profiles');
% K-means Centroid Plot
figure; [cidx, ctrs] = kmeans(expvalues,n,'dist','corr',... %Input
    'rep',5,'disp','final');
for c = 1:16
    subplot(4,4,c);plot(times,ctrs(c,:)'); axis tight,axis off;
end
suptitle('K-Means Centroid Clustering of Profiles');
%% Clustergram %Genes Hierarchical Clustering in Yeast Expression
clustergram(expvalues(:,2:end),'RowLabels',cellstr(genes)','ColumnLabels',times(2:end))
% 1. Input expression 2. genes must be cell array of string
input('Press enter to terminate program'), close all;