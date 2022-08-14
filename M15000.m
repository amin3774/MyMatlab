clc;clear
%tic
%%
PrismCode = fileread("main3.sm");
%main.sm is the prism model created in prism
out1=system('prism MultiAgent.sm -mdp -noprobchecks -exporttrans prism_with_action.txt');%main.sm must be in the current folder of Matlab
if out1  %report an error if prism file is not generated.
   error('"prism_with_action.txt" is not generated.')
end
out1=system('prism main3.sm -exporttransdotstates prism_tuple.txt');
if out1  %report an error if prism file is not generated.
   error('"prism_tuple.txt" is not generated.')
end
clear out1
%%
fileID1 = fopen('prism_with_action.txt','r');
C = textscan(fileID1,'%f %f %f %f %s');
fclose(fileID1);

FlatModelState=C{1,1}(1);%Number of states in the combined model.
FlatModelTran=C{1,3}(1); %Number of transitions in the combined model.
R=sparse(FlatModelState,FlatModelState);%R is the rate matrix
a=C{1,1}(2:end)+1;
b=C{1,3}(2:end)+1;
idx = sub2ind(size(R),a,b);
R(idx)=C{1,4}(2:end);
actions=unique(C{1,5}(2:end));
clear idx a b
%%
f=load('f2');
f=f.f;
% f=ones(FlatModelTran,1);
% f(9)=3;
%f(163)=7;
%f=load('f');f=f.f;

%% Obtain R1 and R2
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
R1=sparse(C1{1,1}(1),C1{1,1}(1));
s1=C1{1,1}(2:end)-min(C1{1,1}(2:end))+1;
t1=C1{1,3}(2:end)-min(C1{1,3}(2:end))+1;
idx1=sub2ind(size(R1),s1,t1);
R1(idx1)=C1{1,4}(2:end);

fid = fopen('prism_with_action2.txt','r');
C2 = textscan(fid,'%f %f %f %f %s');
R2=sparse(C2{1,1}(1),C2{1,1}(1));
s2=C2{1,1}(2:end)-min(C2{1,1}(2:end))+1;
t2=C2{1,3}(2:end)-min(C2{1,3}(2:end))+1;
idx2=sub2ind(size(R2),s2,t2);
R2(idx2)=C2{1,4}(2:end);
%%

%R=full(R);   % Transform sparse matrix to dense matrix.
flat={C{1,1},C{1,3},C{1,4},C{1,5}};
flat{1,1}=flat{1,1}(2:end);flat{1,2}=flat{1,2}(2:end);
flat{1,3}=flat{1,3}(2:end);flat{1,4}=flat{1,4}(2:end);
flat{1,1}(:,2)=flat{1,2};
flat{1,1}(:,3)=flat{1,3};
flat{1,2}=flat{1,4};
flat={flat{1,1},flat{1,2}};

%%
output3=system('prism modul_1.sm -exporttransdotstates prism_tuple1.txt');

fid = fopen('prism_tuple1.txt','r');
u=cell2mat(textscan( fid, '%f [label="%*f%*cn(%f,%f)"];', 'Headerlines', C1{1,3}(1)+2 ));
fclose(fid);
u=u(:,[2 3 1]);
u(:,3)=u(:,3)+1;
%%

fileID2 = fopen('prism_tuple.txt','r');
Tuple=cell2mat(textscan( fileID2, '%f [label="%*f%*cn(%f,%f,%f)"];', 'Headerlines', FlatModelTran+2 ));%Pre-processing. For this Tandem example, last column of "Tuple" is added.
% u=unique(Tuple(:,[2 3]),'row');
% [row,~]=size(u);
 maximum=max(u(:,[1 2]),[],'all');
 m1=sparse(maximum+1,maximum+1);
% u=[u,[1:row]'];
idx=sub2ind(size(m1),u(:,1)+1,u(:,2)+1);
m1(idx)=u(:,3);
idx2=sub2ind(size(m1),Tuple(:,2)+1,Tuple(:,3)+1);
Tuple2=Tuple(:,1);
Tuple2(:,2)=m1(idx2);
Tuple2(:,3)=Tuple(:,4);%Pre-processed 
Tuple=Tuple2;
%fclose(fileID2);
A=Tuple(flat{1,1}(:,1)+1,2:3);
A(:,3)=flat{1,1}(:,3);
A(:,4:5)=Tuple(flat{1,1}(:,2)+1,2:3);
%%
startIndex = regexp(PrismCode,'system');
endIndex = regexp(PrismCode,'endsystem');
sync = PrismCode(startIndex+7:endIndex-1);
startIndex = regexp(sync,'[');
endIndex = regexp(sync,']');
sync=sync(startIndex+1:endIndex-1);
sync=strsplit(sync,{' ',','},'CollapseDelimiters',true);
sync=sync(~cellfun('isempty',sync));

%%
opt=0;
for p=1:200
    i=1;
    f1=f(p,:);
Af=[A,f1']; % f vector is added as the last column of A.
S.a=Af;
S.b=C{1,5}(2:end);%Add action names to Af
S.a(:,7)=S.a(:,6);
S.a(:,8)=[1:FlatModelTran];
if min(S.a(:,2))==0 %%if state numbering in Prsim is started from 0 then we need to change it to start from 1
    S.a(:,2)=S.a(:,2)+1;
    S.a(:,5)=S.a(:,5)+1;
end
if min(S.a(:,1))==0
    S.a(:,1)=S.a(:,1)+1;
    S.a(:,1)=S.a(:,1)+1;
end
indice=find(S.a(:,6)~=1);
amin=indice;
indiceLength=length(indice);

       % R1=[0 2 0 0;0 0 3 4;0 0 0 5;7 0 6 0];
       % R2=[0 0 8 0;9 0 0 10;0 11 0 0 ;12 0 0 0];
        R1orig=R1;R2orig=R2;% Since R1 and R2 will be updated later, I save them  in "R1orig" and "R2orig".

%%

while ~isempty(indice)
    
    if ismember(S.b{indice(i),1},sync) %action "act" is a synchronizing action
        %strfind(S.b,S.b{i,1});   
%         indice2=find(~cellfun(@isempty,strfind(S.b,S.b{indice(i),1})));% find other transitions with same action "act"
%         %%%%%%%%%%%In this part we find the number of "act"-transitions in each
%         %%%%%%%%%%%of the components.
%         W=S.a(indice2,1);
%         W(:,2)=S.a(indice2,4);
%         Y=S.a(indice2,2);
%         Y(:,2)=S.a(indice2,5);
%         W=unique(W,'row'); %Indices of "act"-transitions in c1
%         Y=unique(Y,'row'); %Indices of "act"-transitions in c2
%         [num_action1,col1]=size(W);
%         [num_action2,col2]=size(Y);
%         n=max(num_action1,num_action2); %number of "act" transitions in one of the component(n>=m).
%         m=min(num_action1,num_action2); %number of "act" transitions in the other component.
        %%%%%%%%%%%%%

%         R1=[0 2 0;4 0 3;5 0 0];   %"TEMPORARY"
%         R2=[0 2 0;0 0 3;5 0 0];   %"TEMPORARY"
% R1=[0 1 0 0 0 0 0 0 0 ;0 0 1 0 0 0 0 0 0 ;0 0 0 1 0 0 0 0 0 ;0 0 0 0 1 0 0 0 0 ;0 0 0 0 0 1 0 0 0 ;0 0 0 0 0 0 1 0 0 ;0 0 0 0 0 0 0 1 0 ;0 0 0 0 0 0 0 0 1 ;1 0 0 0 0 0 0 0 0 ];
% R2=[0 1 0 0 0 0 0 0 0 ;0 0 1 0 0 0 0 0 0 ;0 0 0 1 0 0 0 0 0 ;0 0 0 0 1 0 0 0 0 ;0 0 0 0 0 1 0 0 0 ;0 0 0 0 0 0 1 0 0 ;0 0 0 0 0 0 0 1 0 ;0 0 0 0 0 0 0 0 1 ;1 0 0 0 0 0 0 0 0 ];
        matris=S.a(indice2,:);
                  result=gurtest(matris);%%GURTEST%%
                  if strcmp(result.status,'OPTIMAL')
                      opt=opt+1;
%                       uniq1=unique(matris(:,[1 4]),'rows','stable');
%                       idx4=sub2ind(size(R1),uniq1(:,1),uniq1(:,2));
%                       R1(idx4)=result.x(1:length(idx4));
%                       
%                       uniq2=unique(matris(:,[2 5]),'rows','stable');
%                       idx5=sub2ind(size(R2),uniq2(:,1),uniq2(:,2));
%                       R2(idx5)=result.x(length(idx4)+1:end);
%                       disp("R1 and R2 are updated (Gurobi) for SYNC action. ")
                  else
                      newsolution=gurtest7(matris,full(R1orig),full(R2orig));
                      fac1=matris(:,7);
                      fac2=newsolution.x(num_action1+num_action2+1:end);
                      perc(:,p)=100*abs(fac1-fac2)./fac1;
                      mean_perc(p)=mean(perc(:,p));
                      standard(p)=std(perc(:,p));
                      gap(p)=newsolution.mipgap;
                      normf(p)=norm(fac1-fac2);
                      newf(:,p)=fac2;
%                       uniq1=unique(matris(:,[1 4]),'rows','stable');[v1,~]=size(uniq1);
%                       uniq2=unique(matris(:,[2 5]),'rows','stable');[v2,~]=size(uniq2);
%                       y=gurtest4(matris);
%                       %y.x=round(y.x,5,'significant');%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %                       F=nan(FlatModelTran,1);
% %                       F(matris(:,8))=y.x(v1+v2+1:end)
%                       disp(['f[',num2str(matris(:,8)'),']=',num2str(y.x(v1+v2+1:end)')])
%                       disp(num2str(y.x(v1+v2+1:end)'))
%               
%                       error(['Lifting Model Repair info is Impossible(GUROBI) for NON-SYNC Action. \nOne possible modification factor is:\nf('...
%                           ,mat2str(matris(matris(:,7)==1,8)),')= \n',mat2str(y.x(v1+v2+1:end)')],class(FlatModelTran))%add result status to the error message
                  end
               
            S.a(indice2,6)=1;
            indice=find(S.a(:,6)~=1);
           % i=i+1;
      
    else%if "act" is non synchronizing
        
         indice2=find(~cellfun(@isempty,strfind(S.b,S.b{indice(i),1})));% find other transitions with same action "act"
         % prerequisite variables for local repair C1
         ali1=S.a(indice2,[1 4 8]);
         Ind1=find(ismember(ali1(:,[1 2]),S.a(indice(i),[1 4]),'rows'));
         SameSourceInd=ali1(Ind1,3);
         % prerequisite variables for local repair C2
         ali2=S.a(indice2,[2 5 8]);
         Ind2=find(ismember(ali2(:,[1 2]),S.a(indice(i),[2 5]),'rows'));
         SameTargetInd=ali2(Ind2,3);
         % prerequisite variables for Impossibility condition Line58 of
         % Algorithm.
         ali3=S.a(indice2,[1 2 8]);
         %Ind3=find(ismember(ali3(:,[1 2]),S.a(indice(i),[1 2]),'rows'));
         [Ind4,~]=size(unique(ali3(:,[1 2]),'rows'));
         [Ind5,~]=size(ali3);

         if abs(S.a(SameSourceInd,6)-S.a(indice(i),6)*ones(length(SameSourceInd),1))<1e-2 & S.a(indice(i),1)~=S.a(indice(i),4)%local repair in c1 
            R1(S.a(indice(i),1),S.a(indice(i),4))=R1(S.a(indice(i),1),S.a(indice(i),4))*S.a(indice(i),6);%update R1
            S.a(SameSourceInd,6)=1;
            indiceLength=length(indice);
            disp('--> Local Repair in C1 <--')
         elseif abs(S.a(SameTargetInd,6)-S.a(indice(i),6)*ones(length(SameTargetInd),1))<1e-2 & S.a(indice(i),2)~=S.a(indice(i),5)%local repair in c2
            R2(S.a(indice(i),2),S.a(indice(i),5))=R2(S.a(indice(i),2),S.a(indice(i),5))*S.a(indice(i),6);%update R2
            S.a(SameTargetInd,6)=1;
            indiceLength=length(indice);
            disp('--> Local Repair in C2 <--')
            
         elseif Ind4~=Ind5
                arash=S.a(indice2,8);
                ali=unique(S.a(indice2,[1 4]),'rows');
                ali2=ali(ali(:,1)~=ali(:,2),[1 2]); %index of non-sync action in C1
                [row,~]=size(ali2);
                for count=1:row
                    ind(count,:)=find(ismember(S.a(indice2,[1 4]),ali2(count,:),'rows'));
                    Pf1(count,:)=arash(ind(count,:));
                end
                reza=unique(S.a(indice2,[2 5]),'rows');
                ali2=reza(reza(:,1)~=reza(:,2),[1 2]); %index of non-sync action in C1
                [row,~]=size(ali2);
                for count=1:row
                    ind(count,:)=find(ismember(S.a(indice2,[2 5]),ali2(count,:),'rows'));
                    Pf2(count,:)=arash(ind(count,:));
                end
              error(['---> Impossible Condition Satisfied for action "',char(S.b(indice2(1))),...
                  '"(Algorithm Line 58).\n ONLY LOCAL REPAIR IS POSSIBLE! \nPossible local repair solution in C1:\nf('...
                  ,mat2str(Pf1),')=mod_factor \nPossible local repair solution in C2: \nf(',mat2str(Pf2),')=mod_factor<---'],class(row))
         else% In this case we have already checked local repair and impossibility condition and they are not satisfied,
             % so we solve the nonlinear equations directly.
             
                   matris=S.a(indice2,:);
                  result=gurtest(matris);%%GURTEST%%
                  if strcmp(result.status,'OPTIMAL')
                      uniq1=unique(matris(:,[1 4]),'rows','stable');
                      idx4=sub2ind(size(R1),uniq1(:,1),uniq1(:,2));
                      R1(idx4)=result.x(1:length(idx4));
                      
                      uniq2=unique(matris(:,[2 5]),'rows','stable');
                      idx5=sub2ind(size(R2),uniq2(:,1),uniq2(:,2));
                      R2(idx5)=result.x(length(idx4)+1:end);
                      disp("R1 and R2 are updated (Gurobi) for non-sync action. ")
                  else
                      uniq1=unique(matris(:,[1 4]),'rows','stable');[v1,~]=size(uniq1);
                      uniq2=unique(matris(:,[2 5]),'rows','stable');[v2,~]=size(uniq2);
                      y=gurtest4(matris);
%                       F=nan(FlatModelTran,1);
%                       F(matris(:,8))=y.x(v1+v2+1:end)
                      disp(['f[',num2str(matris(:,8)'),']=',num2str(y.x(v1+v2+1:end)')])
                      disp(num2str(y.x(v1+v2+1:end)'))
              
                      error(['Lifting Model Repair info is Impossible(GUROBI). \nOne possible modification factor is:\nf('...
                          ,mat2str(matris(matris(:,7)==1,8)),')= \n',mat2str(y.x(v1+v2+1:end)')],class(FlatModelTran))%add result status to the error message
                  end
           sync=[sync S.b{indice(i)}];       
           S.a(indice2,6)=1;
           indice=find(S.a(:,6)~=1);
         end
          % i=i+1;
    end %if ismember(S.b{indice(i),1},sync)
end %for i=1:length(indice)
end



