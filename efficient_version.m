clc;clear
%main.sm is the prism model created in prism
out1=system('prism main3.sm -mdp -noprobchecks -exporttrans prism_with_action.txt');%main3.sm must be in the current folder of Matlab
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
%R=full(R);   % Transform sparse matrix to dense matrix.
%%
fileID2 = fopen('prism_tuple.txt','r');
Tuple=cell2mat(textscan( fileID2, '%f [label="%*f%*cn(%f,%f)"];', 'Headerlines', FlatModelTran+2 ));
%fclose(fileID2);
%%
flat={C{1,1},C{1,3},C{1,4},C{1,5}};
flat{1,1}=flat{1,1}(2:end);flat{1,2}=flat{1,2}(2:end);
flat{1,3}=flat{1,3}(2:end);flat{1,4}=flat{1,4}(2:end);
flat{1,1}(:,2)=flat{1,2};
flat{1,1}(:,3)=flat{1,3};
flat{1,2}=flat{1,4};
flat={flat{1,1},flat{1,2}};
A=Tuple(flat{1,1}(:,1)+1,2:3);
A(:,3)=flat{1,1}(:,3);
A(:,4:5)=Tuple(flat{1,1}(:,2)+1,2:3);
% %First 2 columns of A are source tuple. Column 3 is the rate.
% %Column 4 and 5 are the destination tuple.
% Astr=string(A);
% Astr(:,6)=flat{1,2};
% %Astr(:,[i j]) = Astr(:,[j i]); %swap columns i, j of a matrix
% 
% var3(1:FlatModelTran,1)="(";
% var4(1:FlatModelTran,1)=",";
% Astr(:,1)=strcat(var3,Astr(:,1),var4);
% var5(1:FlatModelTran,1)=")";
% Astr(:,2)=strcat(Astr(:,2),var5);
% var6(1:FlatModelTran,1)="-";
% var7(1:FlatModelTran,1)=">";
% Astr(:,3)=strcat(var6,Astr(:,6),var4,Astr(:,3),var6,var7);
% Astr(:,4)=strcat(var3,Astr(:,4),var4);
% Astr(:,5)=strcat(Astr(:,5),var5);
% Astr=append(Astr(:,1),Astr(:,2),Astr(:,3),Astr(:,4),Astr(:,5));
% clear var3 var4 var5 var6 var7
%%
% %GUI
% test
% close(test)
% %Here, when I run the gui "test", it generates "guioutput.mat".
% %In 'test_OutputFcn' function, I added 'save guioutput' to save 
% % the output of gui in workspace.
% X = load('guioutput');
% for i=1:FlatModelTran
%     X.handles.popupmenu1.String{i}=Astr(i);
%     X.handles.edit1.Value(i)=1;
% end
% 
% %set(X.handles.popupmenu1,'string','(1,2)->(3,5) \n (2,7)->(3,5)' )
% %test('CALLBACK',h)
%%
% f=ones(FlatModelTran,1);
% f(21)=3;f(26)=7;f(60)=3;f(84)=3;f(108)=3;

f=ones(FlatModelTran,1);
f(19)=4;
f(15)=7;
f([23;27;31;35;40;43;46;49;52;55;60;64;68;72;76;80;85;88;91;94;97;100;105;109;113;117;121;125;130;133;136;139;142;145;150;154;158;162;166;170;175;178;181;184;187;190;195;199;203;207;211;215;220;223;226;229;232;235;240;243;246;249;252;255;258;260;262;264;266;268])=...
[1.75010313905003 1.75010314450272 1.75010315309431...
1.75010369100904 0.39997628895957 0.228557879405469 0.1 0.1 0.1...
0.1 0.39997628895957 0.228557879405469 0.1 0.1 0.1 0.1...
0.399974341402261 0.228556766515578 0.1 0.1 0.1 0.1...
0.399976157135363 0.228557804077351 0.1 0.1 0.1 0.1...
0.399974341402261 0.228556766515578 0.1 0.1 0.1 0.1...
0.39997628895957 0.228557805231538 0.1 0.1 0.1 0.1...
0.399974436165033 0.228556820665733 0.1 0.1 0.1 0.1...
0.39997628895957 0.228557805161179 0.1 0.1 0.1 0.1...
0.399974436165033 0.228556820665733 0.1 0.1 0.1 0.1...
0.39997628895957 0.228557879405469 0.1 0.1 0.1 0.1...
0.39997628895957 0.228557879405469 0.1 0.1 0.1 0.1];
%%
%%%finding synchronizing actions by scanning "main.sm" file,
%%%first we find the variable "PrismCode" which specifies the text between
%%%"system" and "endsytem". Then it finds the action names between "[" and "]".
PrismCode = fileread("main3.sm");
startIndex = regexp(PrismCode,'system');
endIndex = regexp(PrismCode,'endsystem');
sync = PrismCode(startIndex+7:endIndex-1);
startIndex = regexp(sync,'[');
endIndex = regexp(sync,']');
sync=sync(startIndex+1:endIndex-1);
sync=strsplit(sync,{' ',','},'CollapseDelimiters',true);
sync=sync(~cellfun('isempty',sync));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Automatic R1 and R2 generator 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%NOT COMPLETE
%                     equal = regexp(PrismCode,'=');
%                     arrow = regexp(PrismCode,'->');
%                     counter=1;
%                     for i=1:length(equal)
%                        if rem(i,2)~=0
%                            source{counter}=PrismCode(equal(i)+1:arrow(counter)-1);
%                            counter=counter+1;
%                        end
%                     end
%                     source=str2double(source);
%                     amin=strsplit(PrismCode(arrow(counter)+2:end),{' ',','},'CollapseDelimiters',true);
% 
%                     regexp(,'%*s%f%*s:')
% source=PrismCode(startIndex1+1:endIndex1-1);
% source=str2double(source);

% PrismCode=string(PrismCode)
% expression = '=(?<source>\d+)->(?<rate>\d+)';
% tokenNames = regexp(PrismCode,expression,'names');

% regexp(str,'=','match')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
Af=[A,f]; % f vector is added as the last column of A.
S.a=Af;
S.b=C{1,5}(2:end);%Add action names to Af
S.a(:,7)=S.a(:,6);
S.a(:,8)=[1:FlatModelTran];
indice=find(S.a(:,6)~=1);
amin=indice;
indiceLength=length(indice);
%%
%We use this code to extract rate matrices R1 and R2.
%We run each module in PRISM separately and get R1 and R2.
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
i=1;
lookup1=[3.6 12.58];lookup2=[-25 15.48]; %look up table initialisation!
while ~isempty(indice)
%for i=1:indiceLength
    
    if ismember(S.b{indice(i),1},sync) %action "act" is a synchronizing action
        %strfind(S.b,S.b{i,1});   
        indice2=find(~cellfun(@isempty,strfind(S.b,S.b{indice(i),1})));% find other transitions with same action "act"
        %%%%%%%%%%%In this part we find the number of "act"-transitions in each
        %%%%%%%%%%%of the components.
        W=S.a(indice2,1);
        W(:,2)=S.a(indice2,4);
        Y=S.a(indice2,2);
        Y(:,2)=S.a(indice2,5);
        W=unique(W,'row'); %Indices of "act"-transitions in c1
        Y=unique(Y,'row'); %Indices of "act"-transitions in c2
        [num_action1,col1]=size(W);
        [num_action2,col2]=size(Y);
        n=max(num_action1,num_action2); %number of "act" transitions in one of the component(n>=m).
        m=min(num_action1,num_action2); %number of "act" transitions in the other component.
        %%%%%%%%%%%%%

%         R1=[0 2 0;4 0 3;5 0 0];   %"TEMPORARY"
%         R2=[0 2 0;0 0 3;5 0 0];   %"TEMPORARY"
% R1=[0 1 0 0 0 0 0 0 0 ;0 0 1 0 0 0 0 0 0 ;0 0 0 1 0 0 0 0 0 ;0 0 0 0 1 0 0 0 0 ;0 0 0 0 0 1 0 0 0 ;0 0 0 0 0 0 1 0 0 ;0 0 0 0 0 0 0 1 0 ;0 0 0 0 0 0 0 0 1 ;1 0 0 0 0 0 0 0 0 ];
% R2=[0 1 0 0 0 0 0 0 0 ;0 0 1 0 0 0 0 0 0 ;0 0 0 1 0 0 0 0 0 ;0 0 0 0 1 0 0 0 0 ;0 0 0 0 0 1 0 0 0 ;0 0 0 0 0 0 1 0 0 ;0 0 0 0 0 0 0 1 0 ;0 0 0 0 0 0 0 0 1 ;1 0 0 0 0 0 0 0 0 ];
        matris=S.a(indice2,:);
        %%%%%%%%%
%         matris=[1,1,1,2,2,11;% for the case that number of "act" transitions in flat model is less than n*m
%                 1,3,1,2,4,2;
%                 1,5,1,2,6,19;
%                 2,3,1,3,4,1;
%                 2,7,1,3,8,17;
%                 4,1,1,5,2,5;         
%                 4,7,1,5,8,15.4545;
%                 6,9,1,7,10,1;
%                 4 5 1 5 6,8.6364];
        %%%%%%%%%%
        
        if length(indice2)==n*m
            %m=1
            if num_action1==1      %m=1 and its the number of "act" transitions in the first component, c1
                idx = sub2ind(size(R2),S.a(indice2,2),S.a(indice2,5));
                R2(idx)=S.a(indice2,6).*S.a(indice2,3);
                R1(S.a(indice2(1),1),S.a(indice2(1),4))=1;
                disp('--->R1 and R2 updated (m=1 & sync action & Transition number=n*m).<---')
            elseif num_action2==1  %m=1 and its the number of "act" transitions in the second component, c2
                idx = sub2ind(size(R1),S.a(indice2,1),S.a(indice2,4));
                R1(idx)=S.a(indice2,6).*S.a(indice2,3);  
                R2(S.a(indice2(1),2),S.a(indice2(1),5))=1;
                disp('--->R1 and R2 updated (m=1 & sync action & Transition number=n*m).<---')
            % if m>1
            elseif num_action1 >= num_action2 %The case that m is the number of "act"-transitionsa in c2
                      [matris2,~] = sortrows(matris,[1,4,2,5]);%reorder
                      for j=1:m-1
                         for k=1:n
                             matris3(j,k)=matris2(m*(k-1)+1,6)*matris2(m*(k-1)+1,3)/(matris2(m*(k-1)+j+1,6)*matris2(m*(k-1)+j+1,3));%verify condition 1 in the paper
                         end
                      end
                      if abs(matris3-matris3(:,1).*ones(m-1,n))<1e-2
                          R2(matris2(1,2),matris2(1,5))=matris3(1,1);
                          R2(matris2(2,2),matris2(2,5))=1;
                          for t=1:m-2
                            R2(matris2(t+2,2),matris2(t+2,5))=matris3(1,1)/matris3(t+1,1);
                          end
                          for t=0:n-1
                            R1(matris2(t*m+1,1),matris2(t*m+1,4))=matris2(t*m+1,6)*matris2(t*m+1,3)/R2(matris2(t*m+1,2),matris2(t*m+1,5));
                          end
                          disp('--->R1 and R2 updated (m != 1 & sync action & Transition number = n*m).<---')
                      else
                         error('First, Model repair is impossible, because of F value of synchronizing action')% later you
                         %need to add the name of the sync action in the arror message
                      end

            elseif num_action1 < num_action2 %The case that m is the number of "act"-transitionsa in c1
                      [matris2,~] = sortrows(matris,[2,5,1,4]); %reordering
                      for j=1:m-1
                         for k=1:n
                             matris3(j,k)=matris2(m*(k-1)+1,6)*matris2(m*(k-1)+1,3)/(matris2(m*(k-1)+j+1,6)*matris2(m*(k-1)+j+1,3));%verify condition 1 in the paper
                         end
                      end
                      if abs(matris3-matris3(:,1).*ones(m-1,n))<1e-2
                          R1(matris2(1,1),matris2(1,4))=matris3(1,1);
                          R1(matris2(2,1),matris2(2,4))=1;
                          for t=1:m-2
                            R1(matris2(t+2,1),matris2(t+2,4))=matris3(1,1)/matris3(t+1,1);
                          end
                          for t=0:n-1
                            R2(matris2(t*m+1,2),matris2(t*m+1,5))=matris2(t*m+1,6)*matris2(t*m+1,3)/R1(matris2(t*m+1,1),matris2(t*m+1,4));
                          end
                          disp('--->R1 and R2 updated (m != 1 & sync action & Transition number = n*m).<---')
                      else
                         error('Second, Model repair is impossible, because of F value of synchronizing action')% later you
                         %need to add the name of the sync action in the arror message
                      end
            end
            %S.a(indice,6)=1;
        else %Here the number of "act" transitions are less than n*m.
            [a1,b1]=size(unique(matris(:,[1 2]),'rows'));
            [a2,b2]=size(unique(matris(:,[4 5]),'rows'));
            if a1==length(indice2) && a2==length(indice2)
                  %If all "act" transitions in flat model have different 
                  %source and target. e.g. (2,1)->(3,2)   
                  %                        (5,3)->(6,4)
                  %                        (7,6)->(8,7)
                  %This is a simple case and I assign value 1 to x_ij
                  %and f to y_pq.
                 idx3 = sub2ind(size(R1),matris(:,1),matris(:,2));
                 R1(idx3)=matris(:,3);
                 idx4 = sub2ind(size(R2),matris(:,4),matris(:,5));
                 R2(idx4)=matris(:,6);
                  disp('--->R1 and R2 updated (m != 1 & sync action & Transition number < n*m).<---')
            else
                  %If some "act" transitions in flat model have similar  
                  %source and target. e.g. (2,1)->(3,2)   
                  %                        (2,1)->(6,4)
                  %                        (7,6)->(8,7)
                  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
                  lookup1=[3.6 12.58];lookup2=[-25 15.48]; %look up table! 
                  [num,~]=size(matris);
                  [~,idx] = sort(matris(:,6),'descend'); % sort just the 6th column
                matris = matris(idx,:); % sort the whole matrix using the sort indices-->Why I am sorting?I only want
                %to bring the transitions with f~=1, to the first rows of "matris".
                  for j=1:num
                      if ismember(matris(j,[2 5]),lookup2,'rows') && ~ismember(matris(j,[1 4]),lookup1,'rows')
                          R1(matris(j,1),matris(j,4))=matris(j,3)*matris(j,6)/R2(matris(j,2),matris(j,5));
                          lookup1=[lookup1;matris(j,1),matris(j,4)];
                      elseif ~ismember(matris(j,[2 5]),lookup2,'rows') && ismember(matris(j,[1 4]),lookup1,'rows')
                          R2(matris(j,2),matris(j,5))=matris(j,3)*matris(j,6)/R1(matris(j,1),matris(j,4));
                          lookup2=[lookup2;matris(j,2),matris(j,5)];
                      elseif ismember(matris(j,[2 5]),lookup2,'rows') && ismember(matris(j,[1 4]),lookup1,'rows')
                          if abs(R1(matris(j,1),matris(j,4)) * R2(matris(j,2),matris(j,5)) - matris(j,6) *  matris(j,3)) > 1e-3
                            error(['Model repair is not possible. Set f(',num2str(matris(j,8)),')=',num2str(R1(matris(j,1),matris(j,4))...
                                * R2(matris(j,2),matris(j,5))/matris(j,3)),' and try again.'])
                          end
                      elseif ~ismember(matris(j,[2 5]),lookup2,'rows') && ~ismember(matris(j,[1 4]),lookup1,'rows')
                          R1(matris(j,1),matris(j,4))=matris(j,6);
                          R2(matris(j,2),matris(j,5))=matris(j,3);
                          lookup1=[lookup1;matris(j,1),matris(j,4)];
                          lookup2=[lookup2;matris(j,2),matris(j,5)];
                      end
                  end
                  disp('--->R1 and R2 updated (m != 1 & sync action & Transition number < n*m).<---')
            end

        end %if length(indice2)==n*m
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
                  '"(Algorithm Line 58) , Only local repair is possible. \nPossible local repair solution in C1:\nf('...
                  ,mat2str(Pf1),')=mod_factor \nPossible local repair solution in C2: \nf(',mat2str(Pf2),')=mod_factor<---'],class(row))
  
         else% In this case we have already checked local repair and impossibility condition and they are not satisfied,
             % so we solve the nonlinear equations directly.
             
                  matris=S.a(indice2,:);
                  [~,idx] = sort(matris(:,6),'descend'); % sort just the 6th column
                  matris = matris(idx,:); % sort the whole matrix using the sort indices-->Why I am sorting?I only want
                  %to bring the transitions with f~=1, to the first rows of "matris".
                  for j=1:length(indice2)
                      if ismember(matris(j,[2 5]),lookup2,'rows') && ~ismember(matris(j,[1 4]),lookup1,'rows')
                          R1(matris(j,1),matris(j,4))=matris(j,3)*matris(j,6)/R2(matris(j,2),matris(j,5));
                          lookup1=[lookup1;matris(j,1),matris(j,4)];
                      elseif ~ismember(matris(j,[2 5]),lookup2,'rows') && ismember(matris(j,[1 4]),lookup1,'rows')
                          R2(matris(j,2),matris(j,5))=matris(j,3)*matris(j,6)/R1(matris(j,1),matris(j,4));
                          lookup2=[lookup2;matris(j,2),matris(j,5)];
                      elseif ismember(matris(j,[2 5]),lookup2,'rows') && ismember(matris(j,[1 4]),lookup1,'rows')
                          if abs(R1(matris(j,1),matris(j,4)) * R2(matris(j,2),matris(j,5)) - matris(j,6) *  matris(j,3)) > 1e-3
                            error(['Model repair is not possible for non-sync action. Set f(',num2str(matris(j,8)),')=',num2str(R1(matris(j,1),matris(j,4))...
                                * R2(matris(j,2),matris(j,5))/matris(j,3)),' and try again.'])
                          end
                      elseif ~ismember(matris(j,[2 5]),lookup2,'rows') && ~ismember(matris(j,[1 4]),lookup1,'rows')
                          R1(matris(j,1),matris(j,4))=matris(j,6);
                          R2(matris(j,2),matris(j,5))=matris(j,3);
                          lookup1=[lookup1;matris(j,1),matris(j,4)];
                          lookup2=[lookup2;matris(j,2),matris(j,5)];
                      end
                  end 
                disp('--->R1 and R2 updated (Non-sync action).<---')
                disp('sync set is updated.')
                sync=[sync S.b{indice(i)}];

             
         end
           S.a(indice2,6)=1;
           indice=find(S.a(:,6)~=1); 
          % i=i+1;
    end %if ismember(S.b{indice(i),1},sync)
end %for i=1:length(indice)

%At the end of synchronising acion part, you need to make sure that you
%dont check different transition with the same action type. 
% 
% PrismCode = fileread("main.sm");
% startIndex = regexp(PrismCode,'=');
% for i=1:length(startIndex)
%    a{i}=strsplit(PrismCode(startIndex(i)+1),{' ',','},'CollapseDelimiters',true);
% end
% a=a(~cellfun('isempty',a));
sync



