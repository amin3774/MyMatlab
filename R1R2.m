
%We use this code to extract rate matrices R1 and R2. 
%First program reads the main prism file, *.sm.
%Then finds "module" keywords in this file.
%Then creates 2 files: "modul_1.sm" and "modul_2.sm" which represent
%the Prism file of each component.
%Then program runs "modul_1.sm" and "modul_2.sm" in Prism and imports the 
%result files: "prism_with_action1.txt" and "prism_with_action2.txt".
% Matrices R1 and R2 are extracted from these two files.
%
%It seems that this is the best way to extract R1 and R2. The negative
%point is that we need to run Prism 2 more times, which is time consuming.

% I tried to extract R1 and R2 from the main Prism file "*.sm". This method
%is only possible if this file has always a same format. 

PrismCode = fileread("main3.sm"); % read the main prism code
modul=regexp(PrismCode,'module');
part1=PrismCode(1:modul(1)-1);
part2=PrismCode(modul(1):modul(2)+5);
part3=PrismCode(modul(3):modul(4)+5);

part1=string(part1);
part2=string(part2);
part3=string(part3);

modul1=part1 + newline + part2;
modul2=part1+newline+part3;
fid = fopen('modul_1.sm','wt');
fprintf(fid, modul1);
fclose(fid);
fid = fopen('modul_2.sm','wt');
fprintf(fid, modul2);
fclose(fid);
output1=system('prism modul_1.sm -mdp -noprobchecks -exporttrans prism_with_action1.txt');
output2=system('prism modul_2.sm -mdp -noprobchecks -exporttrans prism_with_action2.txt');
fid = fopen('prism_with_action1.txt','r');
C1 = textscan(fid,'%f %f %f %f %s');
fclose(fid);
R1=zeros(C1{1,1}(1),C1{1,1}(1));
s1=C1{1,1}(2:end)-min(C1{1,1}(2:end))+1;
t1=C1{1,3}(2:end)-min(C1{1,3}(2:end))+1;
idx1=sub2ind(size(R1),s1,t1);
R1(idx1)=C1{1,4}(2:end);

fid = fopen('prism_with_action2.txt','r');
C2 = textscan(fid,'%f %f %f %f %s');
R2=zeros(C2{1,1}(1),C2{1,1}(1));
s2=C2{1,1}(2:end)-min(C2{1,1}(2:end))+1;
t2=C2{1,3}(2:end)-min(C2{1,3}(2:end))+1;
idx2=sub2ind(size(R2),s2,t2);
R2(idx2)=C2{1,4}(2:end);
