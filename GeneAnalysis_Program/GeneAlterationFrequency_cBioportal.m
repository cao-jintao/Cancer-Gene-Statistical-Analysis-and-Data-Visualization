%% open file
fid = 'Mutated Genes colorectal carcinoma_cBioportal.txt';
mutated_gene = readtable(fid);mutatedgene = table2cell(mutated_gene);
data = {'genes' 'nmut' 'num' 'perfreq'};
for i = 1:numel(data); 
    s.(data{i}) = mutatedgene(:,i);
    s.(data{i}) = string(s.(data{i}));end
nmut = str2double(s.nmut); num = str2double(s.num);
genes = s.genes; subindex=@(A,loc) A(loc);
for i = 1:numel(s.perfreq);
    perfreq(i) = str2double(subindex(strsplit...
        (s.perfreq(i),'%'),1));end
%%
N = 10; re = 'y';a = 1; b = a+N-1;
while re == 'y',
plot(1:N,perfreq(a:b));xlabel('Genes');ylabel('Percentage');
set(gca,'XTick',1:N,'XTickLabel',genes(a:b));
re = input('continue?(y/n): ','s');
close all,a=a+15;b=a+N-1;end 
%%

