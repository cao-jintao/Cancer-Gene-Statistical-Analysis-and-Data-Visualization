%% Classification
fid = 'TCGACRC_clinical.xlsx';sample = readtable(fid);
cat = input('Input category desired: stage,T,N,M stage: ','s');
subindex = @(a,i) a(i); subjname = string(sample{:,'id'}); %subject names
if cat == 'T';colname = 'tStage'; elseif cat=='N';colname = 'nStage';
elseif cat=='M';colname = 'mStage'; else disp('Category unknown'); end %
strarr = string(sample{:,colname});
strarr(strarr=='NA')= strcat(cat,'9');%mask 'NA' 
for i = 1:numel(strarr)
     str = subindex(strsplit(strarr(i),cat),2);
     stagenum(i) = double(string(subindex(char(str),1)));
end
sets=unique(stagenum); nsets = numel(unique(stagenum));
for i = 1:nsets-1;             %EVAL finding tgrade mask
    eval(sprintf('g%d = find(stagenum==%d);',i,sets(i)));
end
%% CDH17 expression vs Clinicopathological
figure; 
for i = 1:nsets-1
        subplot(2,2,i);        
        eval(sprintf('geneexp = cdh17exp_g%d;',i));
        eval(sprintf('nsample = numel(cdh17exp_g%d);',i));
        hgene = histfit(geneexp,numel(geneexp),'kernel');     %Figure 1;
        maxexp(i) = max(get(hgene(2),'XData'));
        prop=get(gca,'Children');
        set(prop(2),'FaceColor',[0 1-0.2*i 0.1+0.15*i]);
        title(strcat('CRC:CDH17,Class:',colname,',n=',num2str(nsample)));
        legend(sprintf('%s: %d',colname,sets(i)))
        ylabel('No. of Samples');xlabel('Normalized Expression Level');
end
for i = 1:nsets-1;
    peakrank = subindex(unique(maxexp),nsets-i)
    line{i} = sprintf('The peak expression is %.2f in tumor grade %d'...
    ,peakrank,sets(peakrank == maxexp))
end
conclusion = strjoin(string(line),'\n');
figure;annotation('textbox',[0 0 1 1],'string', conclusion)


%
