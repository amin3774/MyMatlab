
%ATTENTION: Do NOT use 'X1','X2',... names as component names, they are
%reserved in the code. Forbiden --->    X1|[action]|X2

%ATTENTION: Never use extra parentheses in system, like: (p1|[a]|p2) or (p1)|[a]|p2.

% This code is feeded with a PRISM code which includes the system, endsystem
% structure and synchronisation(s), like: system= p1|[act_set1]|(p2|[act_set2]|p3).
%The output of this program is a tree of system consisted of some nodes.
%Nodes are the system components (p1,p2,...) and the synchonisation sets,
%like |[a,b]|.
clc;clear;close all
PrismCode = fileread("MultiAgent.sm"); % The PRISM code file main.sm is
                                       % the prism model created in prism.
%%
out1=system('prism MultiAgent.sm -mdp -noprobchecks -exporttrans prism_with_action.txt');%main.sm must be in the current folder of Matlab
if out1  %report an error if prism file is not generated.
   error('"prism_with_action.txt" is not generated.')
end
[out1,~]=system('prism MultiAgent.sm -exporttransdotstates prism_tuple.txt');
if out1  %report an error if prism file is not generated.
   error('"prism_tuple.txt" is not generated.')
end
clear out1

%% find the modules in the "system"
startIndex = regexp(PrismCode,'system');
endIndex = regexp(PrismCode,'endsystem');
system1=PrismCode(startIndex+6:endIndex-1);
system1=strtrim(system1);
main_system=system1;
%
r=regexp(system1,'|[')-1;
l=regexp(system1,']|')+1;
system2=system1;
for i=1:length(r)
    system2(r(i):l(i))=' ';
end
modules=regexpi(system2,'\w+','match');
number_of_modules=length(modules);
%% find the sysnchronisations in the system and rename them Xi
for i=1:length(r)
    syncs{i,1}=system1(r(1):l(1));
    LEN(i)=length(syncs{i,1});
    syncs{i,2}=[' ' 'X' num2str(i) ' '];
    system1(r(1):l(1))='';
    system1=[system1(1:r(1)-1),syncs{i,2},system1(r(1):end)];
     syncs{i,3}=[r(1)];
     syncs{i,4}=[r(1)+length(syncs{i,2})-1];
    r=regexp(system1,'|[')-1;
    l=regexp(system1,']|')+1;
    syncs{i,2}(isspace(syncs{i,2})) = [];
end
alls.a=[modules,syncs{:,2}]';%All the modules and syncs
number_of_nodes=length(alls.a);
for i=1:number_of_nodes
    alls.a{i,2}=i; % assign a number to each node
end
alls.child=alls.a(:,1);
%% Find the main root
[right_par,loc_right_par]=regexp(system1,'(','match');
[left_par,loc_left_par]=regexp(system1,')','match');
if length(loc_right_par)~=length(loc_left_par)
    error('There is missing parenthese(s) in the "SYSTEM" string.')
end
if isempty(right_par) && isempty(left_par) %Exapmple of this case: p1|[a]|p2
    alls.a{1,3}=1;alls.a{2,3}=1;alls.a{3,3}=1;
         for i=1:number_of_nodes-number_of_modules
             if contains(system1,alls.a{i+number_of_modules,1})
                 alls.a{i+number_of_modules,3}=0;
                 root_name=alls.a{i+number_of_modules,1};
                 root_num=i+number_of_modules;
                 break
             end
         end 
%          var1=1:3;var1(root_num)=[];
%          alls.child{root_num,2}=var1(1);alls.child{root_num,3}=var1(2);
elseif loc_right_par(1)~=1 && loc_left_par(end)==length(system1)%Exapmple of this case: p1|[a]|(p2|[b]|p3)
         for i=1:number_of_nodes-number_of_modules
             if contains(system1(1:loc_right_par(1)-1),alls.a{i+number_of_modules,1})
                 alls.a{i+number_of_modules,3}=0;
                 root_name=alls.a{i+number_of_modules,1};
                 root_num=i+number_of_modules;
                 break
             end
         end
         var=system1(1:loc_right_par(1)-1-length(root_name));
         for i=1:number_of_modules
             if contains(var,alls.a{i,1})
                 alls.a{i,3}=root_num;
                 break
             end
         end
         subsys{1,1}=system1(loc_right_par(1)+1:end-1);
         subsys{1,2}=root_num;
elseif loc_right_par(1)==1 && loc_left_par(end)~=length(system1)%Exapmple of this case: (p1|[a]|p2)|[b]|p3
         for i=1:number_of_nodes-number_of_modules
             if contains(system1(loc_left_par(end)-1:end),alls.a{i+number_of_modules,1})
                 alls.a{i+number_of_modules,3}=0;
                 root_name=alls.a{i+number_of_modules,1};
                 root_num=i+number_of_modules;
                 break
             end
         end
         var=system1(loc_left_par+1+length(root_name):end);
         for i=1:number_of_modules
             if contains(var,alls.a{i,1})
                 alls.a{i,3}=root_num;
                 break
             end
         end
         subsys{1,1}=system1(2:loc_left_par(end)-1);
         subsys{1,2}=root_num;
elseif loc_right_par(1)==1 && loc_left_par(end)==length(system1)%Exapmple of this case: (p1|[a]|p2)|[b]|(p3|[c]|p4)
    quaran=[];
    for j=0:length(loc_right_par)-1
        for i=1:length(loc_right_par)
            if loc_left_par(i)>loc_right_par(end-j) && ~ismember(i,quaran)
               a(j+1,:)=[loc_right_par(end-j),loc_left_par(i)];
               quaran=[quaran, i];
               break
            end
        end
    end
    a=sortrows(a);
    k=1;
    for i=1:length(loc_right_par)
            if ~(a(i,1)>a(:,1) & a(i,1)<a(:,2))
                main_par(k,:)=a(i,:);
                k=k+1;
            end
    end
    main_par=sortrows(main_par);
    [b,~]=size(syncs);
    for i=1:b
        if syncs{i,3}>max(main_par(1,:)) && syncs{i,4}>max(main_par(1,:)) &&...
                syncs{i,3}<min(main_par(2,:)) && syncs{i,4}<min(main_par(2,:))
            root.set=syncs{i,2};
            root.pos(1)=syncs{i,3};
            root.pos(2)=syncs{i,4};
        end
    end
    root_name=root.set;
     for i=1:number_of_nodes
         if strcmp(root_name,alls.a{i,1})
             alls.a{i,3}=0;
             root_num=i;
             break
         end
     end 
     subsys={system1(main_par(1,1)+1:main_par(1,2)-1) ; ...%%sub systems
        system1(main_par(2,1)+1:main_par(2,2)-1)};
    subsys{1,2}=root_num;subsys{2,2}=root_num;
end
[~,root.indice]=ismember(root.set,alls.child);
%% This while loop finds the roots in every subsystem and proceeds layer by layer the end of the tree
while ~isempty(subsys)
    sys=subsys{1,1};
    [sub_right_par, sub_loc_right_par]=regexp(sys,'(','match');
    [sub_left_par , sub_loc_left_par] =regexp(sys,')','match');
    
    if isempty(sub_right_par) && isempty(sub_left_par) %Exapmple of this case: p1|[a]|p2
             splitted_sys=split(sys);
             for i=1:number_of_nodes-number_of_modules
                 if strcmp(splitted_sys{2,1},alls.a{i+number_of_modules,1})
                     root_name=alls.a{i+number_of_modules,1};
                     root_num=i+number_of_modules;
                     break
                 end
             end 
             alls.a{root_num,3}=subsys{1,2};
             for i=1:number_of_modules
                 if strcmp(splitted_sys{1,1},alls.a{i,1})
                     alls.a{i,3}=root_num;
                 end
                 if strcmp(splitted_sys{3,1},alls.a{i,1})
                     alls.a{i,3}=root_num;
                 end
             end
             subsys(1,:)=[];
    elseif sub_loc_right_par(1)~=1 && sub_loc_left_par(end)==length(sys)%Exapmple of this case: p1|[a]|(p2|[b]|p3)
            splitted_sys= split(sys(1:sub_loc_right_par(1)-2));
            for i=1:number_of_nodes-number_of_modules
                 if strcmp(splitted_sys{2,1},alls.a{i+number_of_modules,1})
                     alls.a{i+number_of_modules,3}=subsys{1,2};
                     root_name=alls.a{i+number_of_modules,1};
                     root_num=i+number_of_modules;
                     break
                 end
             end
             var=sys(1:sub_loc_right_par(1)-length(root_name));
             splitted_var=split(var);
             for i=1:number_of_modules
                 if strcmp(splitted_var{1,1},alls.a{i,1})
                     alls.a{i,3}=root_num;
                     break
                 end
             end
             subsys(1,:)={sys(sub_loc_right_par(1)+1:end-1),root_num};
    elseif sub_loc_right_par(1)==1 && sub_loc_left_par(end)~=length(sys)%Exapmple of this case: (p1|[a]|p2)|[b]|p3
            splitted_sys=split(sys(sub_loc_left_par(end)+2:end)) ;
            for i=1:number_of_nodes-number_of_modules
                 if strcmp(splitted_sys{1,1},alls.a{i+number_of_modules,1})
                     alls.a{i+number_of_modules,3}=subsys{1,2};
                     root_name=alls.a{i+number_of_modules,1};
                     root_num=i+number_of_modules;
                     break
                 end
             end
             var=sys(sub_loc_left_par+length(root_name):end);
             splitted_var=split(var);
             for i=1:number_of_modules
                 if strcmp(splitted_var{2,1},alls.a{i,1})
                     alls.a{i,3}=root_num;
                     break
                 end
             end
             subsys(1,:)={sys(2:sub_loc_left_par(end)-1),root_num};
    elseif sub_loc_right_par(1)==1 && sub_loc_left_par(end)==length(sys)%Exapmple of this case: (p1|[a]|p2)|[b]|(p3|[c]|p4)
        quaran=[];
        for j=0:length(sub_loc_right_par)-1
            for i=1:length(sub_loc_right_par)
                if sub_loc_left_par(i)>sub_loc_right_par(end-j) && ~ismember(i,quaran)
                   a(j+1,:)=[sub_loc_right_par(end-j),sub_loc_left_par(i)];
                   quaran=[quaran, i];
                   break
                end
            end
        end
        a=sortrows(a);
        k=1;
        for i=1:length(sub_loc_right_par)
                if ~(a(i,1)>a(:,1) & a(i,1)<a(:,2))
                    main_par(k,:)=a(i,:);
                    k=k+1;
                end
        end
        main_par=sortrows(main_par);
        %find main root of subsystem
        aram=sys;
        aram(main_par(2,1):main_par(2,2))=[];
        aram(main_par(1,1):main_par(1,2))=[];
        root_name=aram;
        root_name(isspace(root_name)) = [];
         for i=1:number_of_nodes
             if strcmp(root_name,alls.a{i,1})
                 alls.a{i,3}=subsys{1,2};
                 root_num=i;
             end
         end 
         subsys(1,:)=[];
         subsys(end+1,:)={sys(main_par(1,1)+1:main_par(1,2)-1) , root_num};
         subsys(end+1,:)={sys(main_par(2,1)+1:main_par(2,2)-1) , root_num};
    end 
end

%% Final part to sketch the tree with appropriate labels

for i=1:number_of_nodes
    Final(i)=alls.a{i,3};
end
treeplot(Final)
[x, y] = treelayout (Final);
Final2=cell(1,number_of_modules);
Final2(number_of_modules+1:number_of_nodes)=syncs(:,1);
for i = 1: length (x)
    text (x (i), y (i), [Final2{i},alls.a(i,1),i])
end
title(main_system)
main_system
system1
%% Find the actions and the rate matrix of the flat model: R
% "actions" is a cell array, which the first column shows the action names
% and the rest show the nodes which includes that action.
%It only shows the actions which exist in the flat model.
fileID1 = fopen('prism_with_action.txt','r');
C = textscan(fileID1,'%f %f %f %f %s');
len=length(C{1,1});I=[];     %If deadlock states exist in the flat model,
for i=2:len                  %they are shown by self-loops in the flat model 
   if isempty(C{1,5}{i,1})   %with rate 1, and they have no action name.
        I=[I,i];             %Since they have no action name, we need to remove them.
   end
end
       C{1,1}(I)=[];C{1,2}(I)=[];
       C{1,3}(I)=[];C{1,4}(I)=[];C{1,5}(I)=[];
%%%%
fclose(fileID1);
FlatModelState=C{1,1}(1);%Number of states in the flat model.
FlatModelTran=length(C{1,1})-1;
%FlatModelTran=C{1,3}(1); %Number of transitions in the flat model.
R=sparse(FlatModelState,FlatModelState);%R is the rate matrix
a=C{1,1}(2:end)+1;
b=C{1,3}(2:end)+1;
idx = sub2ind(size(R),a,b);
R(idx)=C{1,4}(2:end);
actions=unique(C{1,5}(2:end));
clear idx a b
%
flat={C{1,1},C{1,3},C{1,4},C{1,5}};
flat{1,1}=flat{1,1}(2:end);flat{1,2}=flat{1,2}(2:end);
flat{1,3}=flat{1,3}(2:end);flat{1,4}=flat{1,4}(2:end);
flat{1,1}(:,2)=flat{1,2};
flat{1,1}(:,3)=flat{1,3};
flat{1,2}=flat{1,4};
flat={flat{1,1},flat{1,2}};
flat{1,1}(:,[1,2])=flat{1,1}(:,[1,2])+1;% In prism, state numbering starts 
                                        %from 0, but I'd like to start from 1.
                         
  %% Obtain Rate matrix of the components
modul=regexp(PrismCode,'module');
part1=string(PrismCode(1:modul(1)-1));
j=1;
for i=1:2:length(modul)
        part2=string(PrismCode(modul(i):modul(i+1)+5));
        mod=part1 + newline + part2;
        fid = fopen('mod.sm','wt');
        fprintf(fid, mod);
        fclose(fid);
        [~,out]=system('prism mod.sm -mdp -noprobchecks -exporttrans prism_with_action_Multi.txt');
        for k=1:length(actions)
            pos=regexp(out,'Modules:');
            line=splitlines(out(pos+8:end));
            line=line{1,1};
            line= line(~isspace(line));
            fid = fopen('prism_with_action_Multi.txt','r');
            CD = textscan(fid,'%f %f %f %f %s');
            fclose(fid);
            Rate.(line).(actions{k,1})=zeros(CD{1,1}(1),CD{1,1}(1));

            s=CD{1,1}(2:end)-min(CD{1,1}(2:end))+1;%source
            t=CD{1,3}(2:end)-min(CD{1,3}(2:end))+1;%target
            s_null=[];t_null=[];
            for l=1:length(s)                       
                if ~strcmp(CD{1,5}(l+1),actions{k,1})
                    s_null=[s_null, s(l)];          
                    t_null=[t_null, t(l)];
                end
            end
            idx_null=sub2ind(size(Rate.(line).(actions{k,1})),s_null,t_null);
            idx=sub2ind(size(Rate.(line).(actions{k,1})),s,t);
            Rate.(line).(actions{k,1})(idx)=CD{1,4}(2:end);
            Rate.(line).(actions{k,1})(idx_null)=0;    
            j=j+1;
        end
end
clear modul part1 part2 mod out pos line CD s t idx fid
%% Obtain Tuple
 fid = fopen('prism_tuple.txt', 'r');      % Open source file.
 for i=1:C{1,3}(1)+2                   % Remove first few lines.
      fgetl(fid) ;  
 end
 buffer = fread(fid, Inf) ;                % Read rest of the file.
 fclose(fid);
 fid = fopen('prism_tuple2.txt', 'w');     % Open destination file.
 fwrite(fid, buffer) ;                     % Save to file.
 fclose(fid) ;
 fid = fopen('prism_tuple2.txt', 'r');     % Extract tuples
 for i=1:FlatModelState
     line=fgetl(fid); 
     [~,f1]=regexp(line,'(');
     [~,f2]=regexp(line,')');
     [~,f3]=regexp(line,'[');
     line2=line(f1+1:f2-1);
     line3=line(1:f3-1);
     Tuple(i,:)=[cell2mat(textscan(line3,'%f')) cell2mat(textscan(line2,...
        '%f','Delimiter',','))'];
 end
 Tuple(:,1)=Tuple(:,1)+1;% In prism, state numbering starts from 0, but 
                         % I'd like to start from 1.
 %% Create Structure A 
 %  A.matrix shows the flat model transitions, half of this matrix shows the
 %  current state and the rest shows the target state.
 %  e.g.: 1 1 1 1 2-->2 3 1 1 2
 A.matrix(:,[1:number_of_modules])=Tuple(flat{1,1}(:,1),[2:end]);
 A.matrix(:,[number_of_modules+1:2*number_of_modules])=Tuple(flat{1,1}(:,2),[2:end]);
 A.action=flat{1,2};
 A.rate=flat{1,1}(:,3);
 %% Create the CHILD
[~,~,var3]=find(unique(cell2mat(alls.a(:,3))));
for i=1:length(var3)
    var4=find(cell2mat(alls.a(:,3))==var3(i));
    alls.child{var3(i),2}=var4(1);
    alls.child{var3(i),3}=var4(2);
end
%% obtain alls the children of the syncs
ch=cell2mat(alls.child(:,[2,3]));
for i=1:number_of_nodes-number_of_modules
    [~,~,dadal]=find(ch(i,:));
    variable=length(dadal);
    j=1;
   while ~isempty(dadal)
       dadal2=cell2mat(alls.child(ch(i,j),[2:end]));
       mehran=length(find(ch(i,:)));
       for k=1:length(dadal2)
           ch(i,mehran+k)=dadal2(k);
       end
       variable=length(find(ch(i,j:end)));
       if variable>1
            [~,~,dadal]=find(ch(i,1+j:end));
       else
           dadal=[];
       end
       j=j+1;
   end
end
clear variable
%%
alls.child(number_of_modules+1:end,2)=syncs(:,1);
[~,b]=size(ch);
alls.child(number_of_modules+1:end,3:b+2)=num2cell(ch);
for i=1:length(actions)
%     regexp(alls.child{i+number_of_modules,2},'\w+','match');
     var4=find(contains(alls.child(number_of_modules+1:end,2),actions(i)))';
     actions(i,[2:length(var4)+1])=num2cell(var4+number_of_modules);
end
scope=actions(:,1);
number_of_actions=length(scope);
for i=1:number_of_actions
    var4=find(contains(alls.child(number_of_modules+1:end,2),actions(i)))'...
        +number_of_modules;
    for j=1:length(var4)
        var5=var4;var5(j)=[];
        if ~ismember(var4(j),cell2mat(alls.child(var5,[3:end])))
            R=scope(i,:);
            scope(i,length(R(~cellfun('isempty',R)))+1)=num2cell(var4(j));
        end
    end
end
[a,b]=size(scope);k=1;
[~,c]=size(alls.child);
for i=1:a
  for j=2:b
     if ~cellfun(@isempty,scope(i,j)) 
        % scope1(k,:)=zeros(1,c);
         scope1.a(k,1)=scope(i,j);
         scope1.a(k,2)=scope(i,1);
         scope1.children(k,:)=alls.child(cell2mat(scope(i,j)),3:end);
         k=k+1;
     end
  end
end
scope=scope1;scope.a
scope.children=cell2mat(scope.children);
[num_of_scopes,~]=size(scope.a);
%%
A.f=ones(1,length(A.rate))'; %A.f: define modification factors manually;
A.f(5)=2;%A.f(10)=3;A.f(26)=5;
A.f([1,8,11,14,17,20,23])=2;
Tmod=A.f(A.f~=1);           %Tmod
Tmod_indice=find(A.f~=1);   %Indice of Tmod
   %% CANNOT
    for i=1:number_of_actions  
       current_act=string(actions(i,1));
       for j=1:number_of_modules
          current_mod=string(alls.child(j,1));
          act.(current_act).(current_mod).CANNOT=[];
          [parent,~]=find(cell2mat(alls.child(number_of_modules+1:end,[3 4]))==j);
          parent=parent+number_of_modules; %Initialisation
          old_parent=j;                    %Initialisation
          cond=true;
          while cond
              salar=char(alls.child(parent,2));
              salar=salar(3:end-2);      % Remove |[ and ]|.
              salar(isspace(salar)) = [];% Remove space from the string.
              salar=split(salar,",");    % Split the action names.
              if isempty(find(strcmp(current_act,salar), 1))% if current_act does Not exist in the parent node...
                temp=cell2mat(alls.child(parent,[3 4]));
                temp=temp(temp~=old_parent);
                  if temp<=number_of_modules
                      act.(current_act).(current_mod).CANNOT=...
                          [act.(current_act).(current_mod).CANNOT;temp];
                  else
                      temp2=nonzeros(cell2mat(alls.child(temp,3:end)));
                      temp2=temp2(temp2<=number_of_modules);
                      act.(current_act).(current_mod).CANNOT=...
                          [act.(current_act).(current_mod).CANNOT;temp2];
                  end
              end
              if parent==root.indice
                  cond=false;
              else
                  old_parent=parent;
                 [parent,~]=find(cell2mat(alls.child(number_of_modules+1:end,[3 4]))==parent);
                 parent=parent+number_of_modules;
              end
          end
       end 
    end
    
   %% MUST
    for i=1:number_of_actions %Vaghti dar root local action nabashad bala miravim, vali agar action bashad ham bayad bala berim ham paiin
       current_act=string(actions(i,1));
       for j=1:number_of_modules
          current_mod=string(alls.child(j,1));
          act.(current_act).(current_mod).MUST=[];
          [parent,~]=find(cell2mat(alls.child(number_of_modules+1:end,[3 4]))==j);
          parent=parent+number_of_modules; %Initialisation
          old_parent=j;                    %Initialisation
          condition=true;
          while condition
              salar=char(alls.child(parent,2));
              salar=salar(3:end-2);      % Remove |[ and ]|.
              salar(isspace(salar)) = [];% Remove space from the string.
              salar=split(salar,",");    % Split the action names.
                  if ~isempty(find(strcmp(current_act,salar), 1))% if current_act exists in the parent node, ham mirim bala ham paiin...
                    temp=cell2mat(alls.child(parent,[3 4]));
                    temp=temp(temp~=old_parent);
                      if temp<=number_of_modules%ok
                          act.(current_act).(current_mod).MUST=...
                              [act.(current_act).(current_mod).MUST;temp];
                      else%be paiin
                          children=nonzeros(cell2mat(alls.child(temp,3:end)))';
                          children=children(children<=number_of_modules);
                          for k=1:length(children)
                              cond=true;
                              [sub_parent,~]=find(cell2mat(alls.child(number_of_modules+1:end,[3 4]))==children(k));
                              sub_parent=sub_parent+number_of_modules;
                              while cond
                                  salar=char(alls.child(sub_parent,2));
                                  salar=salar(3:end-2);      % Remove |[ and ]|.
                                  salar(isspace(salar)) = [];% Remove space from the string.
                                  salar=split(salar,",");    % Split the action names.
                                   if isempty(find(strcmp(current_act,salar), 1))
                                   %if cellfun(@isempty,regexp(current_act,alls.child(sub_parent,2), 'once'))
                                       break
                                   elseif sub_parent==parent
                                       cond=false;
                                       act.(current_act).(current_mod).MUST=...
                                         [act.(current_act).(current_mod).MUST;children(k)];
                                   else
                                       [sub_parent,~]=find(cell2mat(alls.child(number_of_modules+1:end,[3 4]))==sub_parent);
                                       sub_parent=sub_parent+number_of_modules;
                                   end
                              end
                          end
                      end
                      %berim bala
                      if parent==root.indice
                          condition=false;
                      end
                      old_parent=parent;
                      [parent,~]=find(cell2mat(alls.child(number_of_modules+1:end,[3 4]))==parent);
                      parent=parent+number_of_modules;
                  else  %We move to the upper node.
                      old_parent=parent;
                      if parent==root.indice
                            condition = false;
                      end
                      [parent,~]=find(cell2mat(alls.child(number_of_modules+1:end,[3 4]))==parent);
                      parent=parent+number_of_modules;
                      cond1=true;%?
                  end   
          end
        end
    end
   %% May
   %If a component is not niether MUST nor CANNOT, it is MAY.
   for i=1:number_of_actions 
       current_act=string(actions(i,1));
       for j=1:number_of_modules
          current_mod=string(alls.child(j,1));
          act.(current_act).(current_mod).MAY=[];
          u=union(act.(current_act).(current_mod).CANNOT,...
              act.(current_act).(current_mod).MUST);
          act.(current_act).(current_mod).MAY=setdiff(1:number_of_modules,u)';
       end
   end

%%
while ~isempty(Tmod)
    ref_t=Tmod_indice(1);                                                %Take the first transition of Tmod (t_hat)
    indice2=find(~cellfun(@isempty,strfind(A.action,A.action{ref_t,1})));% find other transitions with same action "act"
   % indice2=indice2(indice2~=ref_t);      %Remove ref_t from indice2.
    MS_ref=find(A.matrix(ref_t,1:number_of_modules)~=A.matrix(ref_t,...
        number_of_modules+1:end));                                       %Moving Set of ref_t.    
    % In this for loop, IS of ref_t is obtained.
    % Find out which scope the MS belongs to, and all the seq components in
    % that scope are IS of ref_t.
    
    % To find IS_ref.
    for  i=1:num_of_scopes  
        if ismember(A.action{ref_t,1},scope.a{i,2})
            if ismember(MS_ref,scope.children(i,:))
                var=scope.children(i,:);
                IS_ref=var(var<=number_of_modules);
                IS_ref=find(IS_ref);
                break
            end 
        end
    end
    % Check all transitions with the same action name.
    counter=1;
    for i=1:length(indice2)   
        MS=find(A.matrix(indice2(i),1:number_of_modules)~=A.matrix(indice2(i),...
            number_of_modules+1:end));      %Moving Set of t.
        SS=setdiff(1:number_of_modules,MS); %SS is the complement of MS.
        PS=MS; %Initialise
        
        % To find PS.
        for j=1:length(SS)
            for k=1:length(MS)
                Pi=SS(j);Pj=MS(k);
                if ismember(Pi,act.(char(A.action(ref_t))).(alls.child{Pj,1}).MUST)
                   PS=[PS Pi]; 
                end
                if ismember(Pi,act.(char(A.action(ref_t))).(alls.child{Pj,1}).MAY)...
                        && Rate.(char(modules(Pi))).(char(A.action(ref_t)))...
                        (A.matrix(indice2(i),Pi),A.matrix(indice2(i),Pi))~=0
                   PS=[PS Pi];
               end
            end
        end
        % If PS(t) is not a subset of IS(ref_t),then break and take the
        % next transition of 'indice2'.
        if ~all(ismember(PS, IS_ref))
            break
        else
            PS_t.vector(counter,:)=zeros(1,number_of_modules);
            PS_t.vector(counter,[1:length(PS)])=PS;
            PS_t.indice(counter)=indice2(i);
            
            %create the equation in sym format.
            %ThreeD_mat=zeros(number_of_modules,max(A.matrix,[],'all'),max(A.matrix,[],'all'));
            variable=[];
            for j=1:length(PS)
             variable{j}=['X',num2str(PS(j)),'_',num2str(A.matrix(indice2(i),PS(j))),...
                 '_',num2str(A.matrix(indice2(i),PS(j)+number_of_modules))];
             syms(variable{j})
            end
            e=string(join(variable,'*'));
            eqn(counter)=strcat(e,'==',num2str(A.rate(indice2(i))),'*',num2str(A.f(indice2(i))));
            counter=counter+1;
        end
    end
    eqn=eqn';
    eqns=[];
    for i=1:length(eqn)
        eqns=[eqns;eval(eqn(i))]; 
    end
    solve(eqns)
%     for i=1:length(symvar(eqns))
%         a=symvar(eqns);
%         string(a(i))=optimvar(a(i),'LowerBound',0)
%     end
%     X = optimvar('X',length(symvar(eqns)),1,'LowerBound',0);
%     prob = optimproblem;
%     prob.Objective = norm(X-);
%     options = optimoptions('lsqlin','Display','off');
%     sol = solve(prob,'Options',options)
%     result=solve(eqns)
%     any( structfun(@isempty, result) )

    
    
end


 
 
 