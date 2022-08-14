
%ATTENTION: Do NOT use 'X1','X2',... names as component names, they are
%reserved in the code. Forbiden --->    X1|[action]|X2

%ATTENTION: Never use extra parentheses in system, like: (p1|[a]|p2) or (p1)|[a]|p2.

% This code is feeded with a PRISM code which includes the system, endsystem
% structure and synchronisation(s), like: system= p1|[act_set1]|(p2|[act_set2]|p3).
%The output of this program is a tree consisted of some nodes.
%Nodes are the system components (p1,p2,...) and the synchonisation sets,
%like |[a,b]|.
clc;clear;close all
PrismCode = fileread("Delet3.sm"); % The PRISM code file main.sm is
% the prism model created in prism.
%%
out1=system('prism Delet3.sm -mdp -noprobchecks -exporttrans prism_with_action.txt');%main.sm must be in the current folder of Matlab
if out1  %report an error if prism file is not generated.
    error('"prism_with_action.txt" is not generated.')
end
[out1,~]=system('prism Delet3.sm -exporttransdotstates prism_tuple.txt');
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
r1=regexp(system1,'|[')-1;
l1=regexp(system1,']|')+1;
system2=system1;
for i=1:length(r)
    system2(r(i):l(i))=' ';
end
modules=regexpi(system2,'\w+','match');
number_of_modules=length(modules);clear l1
%% find the synchronisations in the system and rename them Xi
for i=1:length(r)
    syncs{i,1}=system1(r(1):l(1));
    LEN(i)=length(syncs{i,1});
    syncs{i,2}=[' ' 'X' num2str(i) ' '];
    system1(r(1):l(1))='';
    system1=[system1(1:r(1)-1),syncs{i,2},system1(r(1):end)];
    syncs{i,3}=[r(1)];
    syncs{i,4}=[r(1)+length(syncs{i,2})-1];
    %syncs{i,5}=[r(1)+length(syncs{i,1})-1];
    r=regexp(system1,'|[')-1;
    l=regexp(system1,']|')+1;
    syncs{i,2}(isspace(syncs{i,2})) = [];
end
for i=1:length(r1)
    syncs{i,5}=[r1(i)];
    syncs{i,6}=[r1(i)+length(syncs{i,1})-1];
end
alls.a=[modules,syncs{:,2}]';%All the modules and syncs
number_of_nodes=length(alls.a);
for i=1:number_of_nodes
    alls.a{i,2}=i; % assign a number to each node
end
alls.child=alls.a(:,1);clear LEN r l r1
%%
r=regexp(main_system,'[')-1;
l=regexp(main_system,']')+1;
SYS_Parts{1}=main_system(1:r(1)-1);
for i=2:length(l)
    SYS_Parts{i}=main_system(l(i-1)+1:r(i)-1);
end
SYS_Parts{length(l)+1}=main_system(l(length(l))+1:end);clear r l
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
    [~,root.indice]=ismember(root_name,alls.child);
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
    [~,root.indice]=ismember(root_name,alls.child);
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
    [~,root.indice]=ismember(root_name,alls.child);
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
    [~,root.indice]=ismember(root.set,alls.child);
end
%% Arrang paranthesis in matrix "c_mat" in the string "main_system".
[~,right_par_main]=regexp(main_system,'(','match');
[~,left_par_main]=regexp(main_system,')','match');
quaran_right=zeros(1,length(right_par_main));
quaran_left=zeros(1,length(left_par_main));
count=1;
while ~all(quaran_right)
    for j=1:length(right_par_main)
        a_temp=right_par_main;
        a_temp(1:j)=[];
        a_temp=a_temp(find(~quaran_right(length(quaran_right)-length(a_temp)+1:end)));
        b=left_par_main(left_par_main>right_par_main(j));
        b=b(find(~quaran_left(length(quaran_left)-length(b)+1:end)));
        if  ~isempty(a_temp) &&  b(1)<a_temp(1) && ~quaran_right(j) && ~quaran_left(find(left_par_main==b(1)))
            %syncs{j,5}=main_system(right_par_main(j):left_par_main(length(quaran_left)-length(b)+1));
            quaran_right(j)=1;
            quaran_left(find(left_par_main==b(1)))=1;
            c_mat(count,[1 2])=[right_par_main(j), b(1)];
            count=count+1;break
        elseif isempty(a_temp)
            quaran_right(j)=1;
            quaran_left(find(left_par_main==b(1)))=1;
            c_mat(count,[1 2])=[right_par_main(j), b(1)];
            count=count+1;break
        end
    end
end
%Create the processes in every sync set.
for j=1:number_of_nodes-number_of_modules
    count=1;temp=[];
    for k=1:size(c_mat,1)
        if  cell2mat(syncs(j,5))>c_mat(k,1)&& cell2mat(syncs(j,5))<c_mat(k,2) && ...
                cell2mat(syncs(j,6))>c_mat(k,1)&& cell2mat(syncs(j,6))<c_mat(k,2)
            temp(count,[1 2])=[c_mat(k,1),c_mat(k,2)];
            count=count+1;
        elseif j==root.indice-number_of_modules
            syncs{j,7}=main_system;
            count=count+1;
        end
    end
    if ~isempty(temp)
        temp=temp( find(min(abs(temp(:,1)-temp(:,2)))),:);
        syncs{j,7}=main_system(temp(1):temp(2));
    end
end
clear children comb a_temp
%% This while loop finds the roots in every subsystem and proceeds layer by layer till the end of the tree.
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
%title(main_system)
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
R=sparse(FlatModelState,FlatModelState);%R is the rate matrix of the combined model.
a=C{1,1}(2:end)+1;
b=C{1,3}(2:end)+1;

for i=1:length(a)
   Trans_Old(i,:)={num2str(a(i)) num2str(b(i)) char(C{1,5}(i+1))} ;
end
% Trans_Old=cellfun(@num2str,Trans_Old,'un',0);
Trans_Old=categorical(Trans_Old);

I=[];
for i=1:length(a)
    if a(i)==b(i)
        I=[I,i];
    end
end
Trans_Old(I,:)=[];% Tran_Old shows the unique non_selfloop trs.
Trans_Old=unique(Trans_Old,'rows','stable');
FlatNonSelfloopUniqueTrs=length(unique([a,b],'rows'))-length(find(a==b));% Tran_Old shows the the number of unique non_selfloop trs.
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
    Single_modul(j)=string(PrismCode(modul(i):modul(i+1)+5));
    mod=part1 + newline + Single_modul(j);
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
        
       % s=CD{1,1}(2:end)-min(CD{1,1}(2:end))+1;%source
        s=CD{1,1}(2:end)+1;%source
        %t=CD{1,3}(2:end)-min(CD{1,3}(2:end))+1;%target
        t=CD{1,3}(2:end)+1;%target
        s_null=[];t_null=[];rate_null=[];
        for l=1:length(s)
            if strcmp(CD{1,5}(l+1),actions{k,1})
                s_null=[s_null, s(l)];
                t_null=[t_null, t(l)];
                rate_null=[rate_null CD{1,4}(l+1)];
            end
        end
        idx_null=sub2ind(size(Rate.(line).(actions{k,1})),s_null,t_null);
        Rate.(line).(actions{k,1})(idx_null)=rate_null;
    end
    j=j+1;
end
clear modul out pos line CD s t idx fid mod
%% Obtain Tuple
fid = fopen('prism_tuple.txt', 'r');      % Open source file.
prism_tuple = fileread("prism_tuple.txt");
%length(regexp(prism_tuple,'\\n','match'))
number_of_lines=length(regexp(prism_tuple,'\n','match'));%Total number of lines in the prism_tuple.txt file.
for i=1:number_of_lines-FlatModelState-1 % Remove first few lines.
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
%  source state and the rest shows the target state.
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
%% Find all the scopes
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
    actions_struc.(char(actions{i,1}))=cell2mat(actions(i,2:end));
end
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
scope=scope1;
scope.children=cell2mat(scope.children);
[~,a]=size(scope.children);
for i=1:number_of_modules
    for j=1:number_of_actions
        a1=find(strcmp(scope.a(:,2),actions(j,1)));
        a2=cell2mat(scope.a(a1,1)');
        [n1,~]=find(ch==i);
        n1=n1+number_of_modules;
        if isempty(intersect(a2,n1)) && any(Rate.(string(alls.a(i,1))).(string(actions(j,1))),'all')
            scope.a=[scope.a;{i  char(actions(j,1))}] ;
            scope.children=[scope.children;i zeros(1,a-1)];
        end
    end
end
scope.a
[num_of_scopes,~]=size(scope.a);
%%
A.f=ones(1,length(A.rate))'; %A.f: define modification factors manually;
%A.f([40])=2;
% A.f([21 22 23 24 29 30 33 34 423 424 429 432])=2; % For action "a" in Test2
 %A.f([6,8,15,17,33,34,36,37,23,12,25,40,35,41])=[10,20,2,3,10,20,2,3,4,5,30,4,5,30];%action c in "MarkusExample.sm"
%A.f([6    11    16    20    25    29    33    36    41    45    49    52 56    59    62    64])=[ 1 2 3 6 5 10 15 30 4 8 12 24 20 40 60 120];%"PollSystemN=5.sm"
%A.f([8,15])=2;
%A.f([8,15,22,28,35,41,47,52,59,65,71,76,82,87,92,96,103,109,115,120,126,131,136,140,146,151,156,160,165,169,173,176,183,189,195,200,206,211,216,220,226,231,236,240,245,249,253,256,262,267,272,276,281,285,289,292,297,301,305,308,312,315,318,320])=[1,2,7,14,5,10,35,70,1,2,7,14,5,10,35,70,1,2,7,14,5,10,35,70,1,2,7,14,5,10,35,70,1,2,7,14,5,10,35,70,1,2,7,14,5,10,35,70,1,2,7,14,5,10,35,70,1,2,7,14,5,10,35,70];%"PollSystemN=7.sm"
%A.f([8,15,22,28,35,41,47,52,59,65,71,76,82,87,92,96,103,109,115,120,126,131,136,140,146,151,156,160,165,169,173,176,183,189,195,200,206,211,216,220,226,231,236,240,245,249,253,256,262,267,272,276,281,285,289,292,297,301,305,308,312,315,318,320])=[1,3,7,14,5,10,35,70,1,2,7,14,5,10,35,70,1,2,7,14,5,10,35,70,1,2,7,14,5,10,35,70,1,2,7,14,5,10,35,70,1,2,7,14,5,10,35,70,1,2,7,14,5,10,35,70,1,2,7,14,5,10,35,70];%"PollSystemN=7.sm"
%%
%This f is the modification factors for the example in the thesis
% A.f([38
%     39
%     41
%     42
%     43
%     44
%    212
%    213
%    215
%    216
%    217
%    218
%    386
%    387
%    389
%    390
%    391
%    392
%    560
%    561
%    563
%    564
%    565
%    566])=[6
% 6
% 6
% 6
% 6
% 6
% 8
% 8
% 8
% 8
% 8
% 8
% 15
% 15
% 15
% 15
% 15
% 15
% 30
% 30
% 30
% 30
% 30
% 30
% ];
% %%
% A.f([ 11
%     12
%     15
%     21
%     32
%     42
%    170
%    171
%    174
%    180
%    191
%    201
%    329
%    330
%    333
%    339
%    350
%    360
%    488
%    489
%    492
%    498
%    509
%    519])=[6
% 6
% 6
% 6
% 6
% 6
% 8
% 8
% 8
% 8
% 8
% 8
% 15
% 15
% 15
% 15
% 15
% 15
% 30
% 30
% 30
% 30
% 30
% 30
% ];

A.f([1
     3
     4
     8
    12
    14
    16
    17
    21
    25
    55
    57
    58
    61
    65
    67
    69
    70
    73
    77
])=[3
1
4
2
2
4.5
2
5.5
3
3
6
2
8
4
4
9
7
8
6
6
 ];

%%
% A.f([1,2,3,4,18,19,33,34,44,45,46,47,61,62,76,77,87,88,98,110])=6;
% A.f([7,8,9,10,13,14,15,16,24,25,29,30,37,38,41,42,50,51,52,53,56,57,58,59,67,68,72,73,80,81,84,85,91,92,95,96,103,107,113,116])=2; %action "b" in Test2
% A.f([118,119,120,121,132,133,144,145,152,153,154,155,166,167,178,179,186,187,194,203,208,209,210,211,225,226,236,237,247,248,249,250,264,265,275,276,286,287,297,305])=3;
%  rng(1);
%  Index = find(contains(A.action,'loop1a'));
%  A.f(Index)=randi([2,3],size(Index));% Set modification factors randomly, either 2 or 3.
 %% Poll system N=10
% load Satisfiying_f_Poll_10.mat rhs_NEW 
% A.f(Index)=rhs_NEW;
%% Poll system N=11
% Index = find(contains(A.action,'loop1a'));
% % rng(1);
% % A.f(Index)=randi([2,3],size(Index));% Set modification factors randomly, either 2 or 3.
% 
% load Satisfiying_f_Poll_11.mat rhs_NEW 
% A.f(Index)=rhs_NEW;
%% Poll system N=12
% Index = find(contains(A.action,'loop1a'));
% load Satisfiying_f_Poll_12.mat rhs_NEW 
% A.f(Index)=rhs_NEW;
% load ('indN12.mat')
% load ('fN12.mat')
% A.f(indN12)=fN12;
%%
Tmod=A.f(A.f~=1);           %Tmod
% Tmod=1;%remove it
Tmod_indice=find(A.f~=1);   %Indice of Tmod
% Tmod_indice=1;%remove it
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
%If a component is niether MUST nor CANNOT, it is MAY.
for i=1:number_of_actions
    current_act=string(actions(i,1));
    for j=1:number_of_modules
        current_mod=string(alls.child(j,1));
        act.(current_act).(current_mod).MAY=[];
        u=union(act.(current_act).(current_mod).CANNOT,...
            act.(current_act).(current_mod).MUST);
        act.(current_act).(current_mod).MAY=setdiff([1:j-1,j+1:number_of_modules],u)';
    end
end

%%
rem=true(size(A.matrix,1),1);
while ~isempty(Tmod)
    i=1;
    ref_t=Tmod_indice(i);%Take the first transition of Tmod (t_hat).
    if ~rem(ref_t)% if this transition is already considered in previous iterations, skip it.
        continue
    end
    indice2=find(strcmp(A.action,A.action{ref_t,1}));%Find other transitions with same action "act".
    % indice2=indice2(indice2~=ref_t);                %Remove ref_t from indice2.
    MS_ref=find(A.matrix(ref_t,1:number_of_modules)~=A.matrix(ref_t,...
        number_of_modules+1:end));       %Moving Set of ref_t.
    % In the following for loop, IS of ref_t is obtained.
    % Find out the scope which MS belongs to, and all the seq components in
    % that scope are 'IS of ref_t'.
    %% SS_ref and PS_ref
    SS_ref=setdiff(1:number_of_modules,MS_ref);
    PS_ref=MS_ref; %Initialise
    PS_candid_ref=[];
    % To find PS.
    for j=1:length(SS_ref)
        for k=1:length(MS_ref)
            Pi=SS_ref(j);Pj=MS_ref(k);
            if ismember(Pi,act.(char(A.action(ref_t))).(alls.child{Pj,1}).MUST)
                PS_ref=[PS_ref Pi];
            end
            if ismember(Pi,act.(char(A.action(ref_t))).(alls.child{Pj,1}).MAY)...
                    && Rate.(char(modules(Pi))).(char(A.action(ref_t)))(A.matrix(indice2(i),Pi),A.matrix(indice2(i),Pi))~=0
                PS_candid_ref=[PS_candid_ref Pi];
            end
        end
    end
    PS_ref=unique(PS_ref);
    PS_candid_ref=unique(PS_candid_ref);
    cannots_ref=[];
    for j=1:length(MS_ref)
        cannots_ref=[cannots_ref; act.(char(A.action(ref_t))).(alls.child{MS_ref(j),1}).CANNOT];
    end
    PS_candid_ref=setdiff(PS_candid_ref,cannots_ref);
    PS_candid_ref=setdiff(PS_candid_ref,PS_ref);
    %There might be a component which is in a may neighborhood of a
    %moving component AND at the same time in a must neghborhood of another moving component.
    %This component is for sure participating and that is why it is removed from PS_candid.
    %% Check PS_candid components
    for j=1:length(PS_candid_ref)
        % Find the smallest subtree,"The_node" consisting of PS_candid and a moving
        % component
        Bol=true;
        %"Bol" is a boolean varaible. For every component in PS_candid,
        %if 'Bol' is true at end, then it is a participating component.
        The_node=[];The_node_len=[];
        a1=find(strcmp(actions(:,1),A.action(ref_t)));
        for k=1:number_of_nodes-number_of_modules
            if ismember(number_of_modules+k,cell2mat(actions(a1,[2:end]))) %Check the nodes which incldes ref_action.
                am=cell2mat(alls.child(k+number_of_modules,[3:end]));
                am=am(find(am));
                if ismember(PS_candid_ref(j),am) && ~isempty(intersect(MS_ref,am))
                    The_node=[The_node,k+number_of_modules];
                    The_node_len=[The_node_len,length(find(cell2mat(alls.child(k+number_of_modules,[3:end]))))];
                end
            end
        end
        [~,The_node_len]=min(The_node_len);
        The_node=The_node(The_node_len);
        [parent,~]=find(cell2mat(alls.child(number_of_modules+1:end,[3,4]))==PS_candid_ref(j));
        parent=parent+number_of_modules;
        while parent ~= The_node
            if ismember(parent,cell2mat(actions(a1,[2:end])))
                children=cell2mat(alls.child(parent,[3:end]));
                children=children(find(children));
                children=children(children<=number_of_modules);
                comb=part1 + newline;
                for k=1:length(children)
                    comb=comb + newline + Single_modul(children(k));
                end
                comb=comb + newline + 'system' + newline + syncs{parent-number_of_modules,7} + newline + 'endsystem';
                fid = fopen('mod.sm','wt');
                fprintf(fid, comb);
                fclose(fid);
                system('prism mod.sm -mdp -noprobchecks -exporttrans mod.txt');
                fid = fopen('mod.txt','r');
                C1 = textscan(fid,'%f %*f %f %*f %s');%Extract only source and target states and the action name.
                FlatModelState2=C1{1,1}(1);
                fclose(fid);
                C2.transitions=[C1{1,1}(2:end),C1{1,2}(2:end)];
                C2.action=C1{1,3}(2:end);
                % Obtain Tuple
                system('prism mod.sm -exporttransdotstates mod_tuple.txt');
                fid = fopen('mod_tuple.txt','r');
                prism_tuple = fileread("mod_tuple.txt");
                number_of_lines=length(regexp(prism_tuple,'\n','match'));%Total number of lines in the prism_tuple.txt file.
                for k=1:number_of_lines-FlatModelState2-1 % Remove first few lines.
                    fgetl(fid) ;
                end
                buffer = fread(fid, Inf) ;             % Read rest of the file.
                fclose(fid);
                fid = fopen('mod_tuple.txt', 'w');     % Open destination file.
                fwrite(fid, buffer) ;                  % Save to file.
                fclose(fid) ;
                fid = fopen('mod_tuple.txt', 'r');     % Extract tuples
                for k=1:FlatModelState2
                    line=fgetl(fid);
                    [~,f1]=regexp(line,'(');
                    [~,f2]=regexp(line,')');
                    [~,f3]=regexp(line,'[');
                    line2=line(f1+1:f2-1);
                    line3=line(1:f3-1);
                    Tuple2(k,:)=[cell2mat(textscan(line3,'%f')) cell2mat(textscan(line2,...
                        '%f','Delimiter',','))'];
                end
                Tuple2(:,1)=Tuple2(:,1)+1;% In prism, state numbering starts from 0, but
                % I'd like to start from 1.
                ind=find(ismember( Tuple2(:,2:end),A.matrix(indice2(i),children),'rows'));
                if ismember([ind ind],C2.transitions,'rows')
                    ind2=find(ismember(C2.transitions,[ind ind],'rows'));
                    if any(strcmp(C2.action(ind2),A.action(ref_t)))
                        Bol=true;
                    else
                        Bol=false;
                    end
                else
                    Bol=false;
                end
            end
            [parent,~]=find(cell2mat(alls.child(number_of_modules+1:end,[3,4]))==parent);% The new parent is the parent of the old one.
            parent=parent+number_of_modules;
        end
        if Bol
            PS_ref=[PS_ref, PS_candid_ref(j)];
        end
    end
    %% To find IS_ref.
    IS_ref=PS_ref; %Initialise
    for j=1:length(IS_ref)
        IS_ref=[IS_ref,act.(A.action{ref_t,1}).(alls.child{IS_ref(j),1}).MAY'];
        IS_ref=unique(IS_ref,'stable');
    end
    
    scope_ref=[];
    for  j=1:num_of_scopes
        if ismember(A.action{ref_t},scope.a{j,2})
            if ismember(MS_ref,scope.children(j,:))
                %                 var=scope.children(i,:);
                %                 IS_ref=var(var<=number_of_modules);
                %                 IS_ref=find(IS_ref);
                scope_ref=j;
                break
            end
        end
    end
    %%
    % Check all transitions with the same action name.
    counter=1;eqns=[];
    for c1=1:length(indice2)
        counter2=1;
        MS=find(A.matrix(indice2(c1),1:number_of_modules)~=A.matrix(indice2(c1),...
            number_of_modules+1:end));      %Moving Set of t.
        if isempty(MS) % If there is a self-loop in the combined model, we ignore it.
            continue
        end
        SS=setdiff(1:number_of_modules,MS); %SS is the complement of MS.
        %%
        PS=MS; %Initialise
        PS_candid=[];
        % To find PS.
        for j=1:length(SS)
            for k=1:length(MS)
                Pi=SS(j);Pj=MS(k);
                if ismember(Pi,act.(char(A.action(ref_t))).(alls.child{Pj,1}).MUST)
                    PS=[PS Pi];
                end
                if ismember(Pi,act.(char(A.action(ref_t))).(alls.child{Pj,1}).MAY)...
                        && Rate.(char(modules(Pi))).(char(A.action(ref_t)))(A.matrix(indice2(c1),Pi),A.matrix(indice2(c1),Pi))~=0
                    PS_candid=[PS_candid Pi];
                end
            end
        end
        PS=unique(PS);
        PS_candid=unique(PS_candid);
        cannots=[];
        for j=1:length(MS)
            cannots=[cannots; act.(char(A.action(ref_t))).(alls.child{MS(j),1}).CANNOT];
        end
        PS_candid=setdiff(PS_candid,cannots);
        PS_candid=setdiff(PS_candid,PS);
        %There might be a component which is in a may neighborhood of a
        %moving component AND at the same time in a must neghborhood of another moving component.
        %This component is for sure participating and that is why it is removed from PS_candid.
        %% Check PS_candid components
        for j=1:length(PS_candid)
            % Find the smallest subtree,"The_node" consisting of PS_candid and a moving
            % component
            Bol=true;
            %"Bol" is a boolean varaible. For every component in PS_candid,
            %if 'Bol' is true at end, then it is a participating component.
            The_node=[];The_node_len=[];
            a1=find(strcmp(actions(:,1),A.action(ref_t)));
            for k=1:number_of_nodes-number_of_modules
                if ismember(number_of_modules+k,cell2mat(actions(a1,[2:end]))) %Check the nodes which incldes ref_action.
                    am=cell2mat(alls.child(k+number_of_modules,[3:end]));
                    am=am(find(am));
                    if ismember(PS_candid(j),am) && ~isempty(intersect(MS,am))
                        The_node=[The_node,k+number_of_modules];
                        The_node_len=[The_node_len,length(find(cell2mat(alls.child(k+number_of_modules,[3:end]))))];
                    end
                end
            end
            [~,The_node_len]=min(The_node_len);
            The_node=The_node(The_node_len);
            [parent,~]=find(cell2mat(alls.child(number_of_modules+1:end,[3,4]))==PS_candid(j));
            parent=parent+number_of_modules;
            while parent ~= The_node
                if ismember(parent,cell2mat(actions(a1,[2:end])))
                    children=cell2mat(alls.child(parent,[3:end]));
                    children=children(find(children));
                    children=children(children<=number_of_modules);
                    comb=part1 + newline;
                    for k=1:length(children)
                        comb=comb + newline + Single_modul(children(k));
                    end
                    comb=comb + newline + 'system' + newline + syncs{parent-number_of_modules,7} + newline + 'endsystem';
                    fid = fopen('mod.sm','wt');
                    fprintf(fid, comb);
                    fclose(fid);
                    system('prism mod.sm -mdp -noprobchecks -exporttrans mod.txt');
                    fid = fopen('mod.txt','r');
                    C1 = textscan(fid,'%f %*f %f %*f %s');%Extract only source and target states and the action name.
                    FlatModelState2=C1{1,1}(1);
                    fclose(fid);
                    C2.transitions=[C1{1,1}(2:end),C1{1,2}(2:end)]+1;% In prism, state numbering starts from 0, but
                    % I'd like to start from 1. (The reason for +1)
                    C2.action=C1{1,3}(2:end);
                    % Obtain Tuple
                    system('prism mod.sm -exporttransdotstates mod_tuple.txt');
                    fid = fopen('mod_tuple.txt','r');
                    prism_tuple = fileread("mod_tuple.txt");
                    number_of_lines=length(regexp(prism_tuple,'\n','match'));%Total number of lines in the prism_tuple.txt file.
                    for k=1:number_of_lines-FlatModelState2-1 % Remove first few lines.
                        fgetl(fid) ;
                    end
                    buffer = fread(fid, Inf) ;             % Read rest of the file.
                    fclose(fid);
                    fid = fopen('mod_tuple.txt', 'w');     % Open destination file.
                    fwrite(fid, buffer) ;                  % Save to file.
                    fclose(fid) ;
                    fid = fopen('mod_tuple.txt', 'r');     % Extract tuples
                    for k=1:FlatModelState2
                        line=fgetl(fid);
                        [~,f1]=regexp(line,'(');
                        [~,f2]=regexp(line,')');
                        [~,f3]=regexp(line,'[');
                        line2=line(f1+1:f2-1);
                        line3=line(1:f3-1);
                        Tuple2(k,:)=[cell2mat(textscan(line3,'%f')) cell2mat(textscan(line2,...
                            '%f','Delimiter',','))'];
                    end
                    Tuple2(:,1)=Tuple2(:,1)+1;% In prism, state numbering starts from 0, but
                    % I'd like to start from 1.
                    ind=find(ismember( Tuple2(:,2:end),A.matrix(indice2(c1),children),'rows'));
                    if ismember([ind ind],C2.transitions,'rows')
                        ind2=find(ismember(C2.transitions,[ind ind],'rows'));
                        if any(strcmp(C2.action(ind2),A.action(ref_t)))
                            Bol=true;
                        else
                            Bol=false;
                            c1
                            break
                        end
                    else
                        Bol=false;
                        c1
                        break
                    end
                end
                [parent,~]=find(cell2mat(alls.child(number_of_modules+1:end,[3,4]))==parent);% The new parent is the parent of the old one.
                parent=parent+number_of_modules;
            end
            if Bol
                PS=[PS, PS_candid(j)];
            end
        end
        
        A.PS{indice2(c1)}=PS';
        %%
        % If PS(t) is not a subset of IS(ref_t),then skip this t and take the
        % next transition of 'indice2'.
        if ~all(ismember(PS, IS_ref))
            continue
        end
        rem(indice2(c1))=false;
        PS_t.vector(counter,:)=zeros(1,number_of_modules);
        PS_t.vector(counter,[1:length(PS)])=PS;
        PS_t.indice(counter)=indice2(c1);
        
        %Obtain the 'variable', where the first column shows the indice
        %of the variable, the second one is the indice
        %of the component, third one shows the source state and the
        %fourth one shows the target state in that component.
        
        if ismember(PS,MS) % if PS=MS
            eqn=[];
            for j=1:length(PS)
                variable(counter2,:)=[PS(j),A.matrix(indice2(c1),PS(j)),...
                    A.matrix(indice2(c1),PS(j)+number_of_modules)];
                eqn=[eqn, variable(counter2,:)];
                vari{j}=['X',num2str(PS(j)),'_',num2str(A.matrix(indice2(c1),PS(j))),...
                    '_',num2str(A.matrix(indice2(c1),PS(j)+number_of_modules))];
                syms(vari{j})
                counter2=counter2+1;
            end
            e=string(join(vari,'*'));
            Eq_string(counter)=strcat(e,'==',num2str(A.rate(indice2(c1))),'*',num2str(A.f(indice2(c1))));
            %                 eqns.lhs(counter,1:length(eqn))=eqn;
            %                 eqns.rhs(counter,:)=A.rate(indice2(i))*A.f(indice2(i));
            counter=counter+1;
        else %(RSLC)If there are some self-loops and |PS|>|MS|.
            %% RSLC: find the smallest subtree consists of PS: variable name:"local_root"
            a=cell2mat(alls.child(number_of_modules+1:end,3:end));
            b=[];c=[];
            for j=1:number_of_nodes-number_of_modules
                if ismember(PS,a(j,:))
                    b=[b, number_of_modules+j];
                    c=[c,length(find(a(j,:)))];
                end
            end
            [~,r]=min(c);
            local_root=b(r);%smallest subtree of PS
            
            b=a(local_root-number_of_modules,:);
            ch=b(find(b));ch=[ch,local_root];
            len=length(ch);
            remain=ch;
            %ch(2,:)=0;
            comb=cell(1,len);
            while ~isempty(remain)
                for j=1:len % check the leafs.
                    if ch(1,j)<=number_of_modules
                        musts=[];
                        for k=1:length(MS)
                            musts=[musts; act.(char(A.action(ref_t))).(alls.child{MS(k),1}).MUST];
                            musts=unique(musts);
                        end
                        if ismember(ch(1,j),PS) && ~ismember(ch(1,j),MS) && ~ismember(ch(1,j),musts)
                            comb{j}=ch(1,j);
                            remain(remain==ch(1,j))=[];
                        else
                            %comb{j}=0;
                            remain(remain==ch(1,j))=[];
                        end
                    end
                end
                for j=1:len % check the nodes.
                    if ch(1,j)>number_of_modules && ~ismember(alls.child{ch(j),3},remain) &&...
                            ~ismember(alls.child{ch(j),4},remain) && ismember(ch(j),remain)
                        left=find(ch==alls.child{ch(j),3});
                        right=find(ch==alls.child{ch(j),4});
                        if ismember(ch(1,j),actions_struc.(char(A.action(ref_t))))
                            if isempty(comb{left}) && isempty(comb{right})
                                %comb{j}=[];
                                remain(remain==ch(1,j))=[];
                            elseif ~isempty(comb{left}) && isempty(comb{right})
                                comb{j}=comb{left};
                                remain(remain==ch(1,j))=[];
                            elseif isempty(comb{left}) && ~isempty(comb{right})
                                comb{j}=comb{right};
                                remain(remain==ch(1,j))=[];
                            elseif ~isempty(comb{left}) && ~isempty(comb{right})
                                mat1=comb{left};
                                mat2=comb{right};
                                ma=size(mat1,1);
                                mb=size(mat2,1);
                                [b1,b2]=ndgrid(1:ma,1:mb);
                                comb{j}=[mat1(b1,:),mat2(b2,:)];
                                remain(remain==ch(1,j))=[];
                            end
                        else %Union of left and right leaves.
                            sx = size(comb{left});
                            sy = size(comb{right});
                            am = max(sx(2),sy(2));
                            comb{j} = [[comb{left},zeros(abs([0 am]-sx))];[comb{right},zeros(abs([0,am]-sy))]];
                            remain(remain==ch(1,j))=[];
                        end
                    end
                end
            end
            % Here we create the equations if |PS|>|MS|, and subsequently there exists
            %rslc set. First, the vaiables corresponding MS, then
            %variables corresponding Q and finally rslc.
            rslc=comb{find(ch==local_root)} ; %RSLC result.
            multi_trans_idx=indice2(find(ismember(A.matrix(indice2,:),A.matrix(indice2(c1),:),'rows')));
            if ~isempty(rslc) && size(rslc,1)~=length(multi_trans_idx)
                error('Error 1 in RSLC function. The number of self-loop combinations is different with the number of multi transitions.')
            end
            if isempty(rslc) && length(multi_trans_idx)>1
                error('Error 2 in RSLC function. RSLC is empty but there are some multi transitions.')
            end
            eqn=[];vari=[];
            for j=1:length(MS) %counter2 is the number of variables.
                variable(counter2,:)=[MS(j),A.matrix(indice2(c1),MS(j)),...
                    A.matrix(indice2(c1),MS(j)+number_of_modules)];
                eqn=[eqn, variable(counter2,:)];
                vari{j}=['X',num2str(MS(j)),'_',num2str(A.matrix(indice2(c1),MS(j))),...
                    '_',num2str(A.matrix(indice2(c1),MS(j)+number_of_modules))];
                syms(vari{j})
                counter2=counter2+1;
            end
            % Q is the set of components defined in line 48 of the algorithm.
            Q=intersect(musts,setdiff(PS,MS));
            for j=1:length(Q)
                variable(counter2,:)=[Q(j),A.matrix(indice2(c1),Q(j)),...
                    A.matrix(indice2(c1),Q(j)+number_of_modules)];
                eqn=[eqn, variable(counter2,:)];
                vari{counter2}=['X',num2str(Q(j)),'_',num2str(A.matrix(indice2(c1),Q(j))),...
                    '_',num2str(A.matrix(indice2(c1),Q(j)+number_of_modules))];
                syms(vari{counter2})
                counter2=counter2+1;
            end
            %rslc
            if ~isempty(rslc) %if rslc is not empty.
                part3={};e={}; part2={};
                for j=1:size(rslc,1)
                    s1=rslc(j,:);
                    s1=s1(find(s1)); %remove zero elements
                    s2=length(s1);
                    for k=1:s2
                        part2{j,k}=['X',num2str(rslc(j,k)),'_',num2str(A.matrix(indice2(c1),rslc(j,k))),...
                            '_',num2str(A.matrix(indice2(c1),rslc(j,k)+number_of_modules))];
                        syms(part2{j,k})
                        counter2=counter2+1;
                    end
                    part3(j,1:length([vari part2(j,:)]))=[vari part2(j,:)];
                    e(j)=join(part3(j,:),'*');
                end
                complete=string(join(e,'+'));
                RHS=0;
                for j=1:size(rslc,1) %Right Hand Side of the equation, which is the product of f and rates.
                    RHS=RHS+A.rate(multi_trans_idx(j))*A.f(multi_trans_idx(j));
                end
                Eq_string(counter)=strcat(complete,'==',num2str(RHS));
                %                 eqns.lhs(counter,1:length(eqn))=eqn;
                %                 eqns.rhs(counter,:)=A.rate(indice2(i))*A.f(indice2(i));
            else % if rslc is empty.
                e=string(join(vari,'*'));
                Eq_string(counter)=strcat(e,'==',num2str(A.rate(indice2(c1))),'*',num2str(A.f(indice2(c1))));
            end
            rem(multi_trans_idx)=false;
            counter=counter+1;
        end
        %end
    end
    Eq_string=Eq_string';
    Eq_string=unique(Eq_string,'stable');
    equations=[];
    for j=1:length(Eq_string)
        equations=[equations;eval(Eq_string(j))];
    end
    assume(symvar(equations)>0);% assume that all the variables are positive and non-zero.
    solution1=solve(equations,symvar(equations), 'IgnoreAnalyticConstraints',1) %Solve system of eqs.
    solution1_num=[];
    if isstruct (solution1)
        solution1_num = struct2cell(solution1);
        solution1_num = cat(2,solution1_num{:});
        solution1_num=solution1_num';
    end
    %% If there exists a solution, then we try to find an optimised solution.
    if ~isempty(solution1_num)
        A.f(~find(rem))=1;
        %Tmod(indice2)=[];
        disp(['A solution is found for action ',A.action{ref_t}])
        break
        %continue
    end
    %         variable=unique(variable,'row');
    %         [num_var,~]=size(variable);
    %         ThreeD_mat=zeros(number_of_modules,max(A.matrix,[],'all'),max(A.matrix,[],'all'));
    %         %    idx=sub2ind(size(ThreeD_mat),variable(:,1),variable(:,2),variable(:,3));
    %         for i=1:num_var
    %             ThreeD_mat(variable(i,1),variable(i,2),variable(i,3))=i;
    %         end
    %         a=1:num_var;a=a';
    %         variable=[a,variable];
    %         [~,a]=size(eqns.lhs);
    %         equation=[];
    %         for i=1:length(eqns.rhs)
    %             b=[];
    %             for j=1:3:a
    %                 if eqns.lhs(i,j)
    %                     b=[b,ThreeD_mat(eqns.lhs(i,j),eqns.lhs(i,j+1),eqns.lhs(i,j+2))];
    %                 end
    %             end
    %             equation(i,1:length(b))=b;
    %         end
    %         eqns.lhs=equation;
    %         X = optimvar('X',num_var,1,'LowerBound',0);
    %         prob = optimproblem;
    %         for i=1:num_var
    %             Y(i)= Rate.(char(modules(variable(i,2)))).(A.action{ref_t,1})(variable(i,3),variable(i,4));
    %         end
    %         prob.Objective = sum((X-Y').^2);
    %         for i=1:length(eqns.rhs)
    %             a=eqns.lhs(i,:);
    %             a=a(find(a)); %Do no remove 'find'!
    %             part=[];
    %             for j=1:length(a)
    %                 part=[part,'X(',num2str(a(j)),')'];
    %                 if j~=length(a)
    %                     part=[part,'*'];
    %                 end
    %             end
    %             part=['prob.Constraints.cons',num2str(i),'=',part,'==',num2str(eqns.rhs(i))];
    %             eval(part);
    %         end
    %         %     options = optimoptions('lsqlin','Display','off');
    %         x0.X=Y;
    %         solution2 = solve(prob,x0);
    %         solution2 = struct2cell(solution2);
    %         solution2 = cat(2,solution2{:});
    %         % sol = solve(prob,x0,'Options',options);
    %%     any( structfun(@isempty, result) )
    % if solution1 is empty, we try to add the action to the other
    % sync sets.
    %'NoActSets' is a vector includes all the sync sets with no
    % ref_action.
    for j=1:number_of_actions
        if strcmp(A.action{ref_t},actions{j,1})
            NoActSets=cell2mat(actions(j,2:end));
            NoActSets=setdiff(number_of_modules+1:number_of_nodes,NoActSets);
            break
        end
    end
    % var2 is a set of the sync sets in side the scope.
    var2=scope.children(scope_ref,:);
    var2=var2(find(var2)); %#ok<*FNDSB>
    var2=var2(var2>number_of_modules);
    var2=intersect(var2,NoActSets);
    % var7 is a set of the sync sets outside of the scope until the main root.
    n=cell2mat(scope.a(scope_ref,1));
    var7=[];
    while n~=root.indice
        [~,a1]=ismember (n,cell2mat(alls.child(number_of_modules+1:end,[3 4])));
        if a1>number_of_nodes-number_of_modules
            a1=a1-(number_of_nodes-number_of_modules);
        end
        n=a1+number_of_modules;
        if ismember(n,NoActSets)
            var7=[var7 n]; % Now var7 vector is ordered, i.e. each element of var7 is the child of the next element.
        end
    end
    %%
    add_action=[];
    Single_modul_PRIME=Single_modul;

    %var2=[25:-1:15]; % N=12
    %var2=[23:-1:14]; % N=11
    %var2=[21:-1:13];  % N=10
    %var2=[19:-1:12]; % N=9
    %var2=[17:-1:11];% N=8
    %var2=[15:-1:10]; % N=7
    %var2=[13:-1:9];  % N=6 
   % var2=[11 13 10];%N=5
   %var2=7;
   for j=1:length(var2)
        old_length_add_action=length(add_action);
        flag2=true;
        %first add the action to one of the sync sets inside the
        %scope. If it does not work or if var2 is empty, try var7,
        %i.e. the sync sets out of the scope.
        current_node=var2(j)-number_of_modules;
        im_children=cell2mat(alls.child(var2(j),[3 4]));%Two immediate children of 'var2(j)'
        B1=im_children(1)<=number_of_modules;
        B2=im_children(2)<=number_of_modules;
        flag=true; % if flag remains true at the end of next "if", then there is no type A sp. tr.
        if B1 && B2%Both immediate children are leaves.
            for k=1:length(Tuple) %Tuple is the list of all reachable states in the flat model.
                first=Tuple(k,im_children(1)+1);
                B3=any(Rate.(char(modules(im_children(1)))).(char(A.action(ref_t)))(first, :));
                second=Tuple(k,im_children(2)+1);
                B4=any(Rate.(char(modules(im_children(2)))).(char(A.action(ref_t)))(second,:));
                if B3 && B4
                    flag=false;%type A sp.tr.
                    break
                end
            end
%             if flag
%                 add_action=[add_action, var2(j)];
%             end
        elseif xor(B1,B2)
            if B1
                leaf=1;
                not_leaf=2;
            else
                leaf=2;
                not_leaf=1;
            end
            %X1=any(Rate.(char(modules(im_children(leaf)))).(char(A.action(ref_t)))(:));%leaf
            %X2:non_leaf; It is a process.
            children=cell2mat(alls.child(im_children(not_leaf),[3:end]));
            children=children(find(children));
            children=children(children<=number_of_modules);
            comb=part1 + newline;
            for k=1:length(children)
                comb=comb + newline + Single_modul_PRIME(children(k));
            end
            comb=comb + newline + 'system' + newline + syncs{im_children(not_leaf)-number_of_modules,7} + newline + 'endsystem';
            fid = fopen('mod.sm','wt');
            fprintf(fid, comb);
            fclose(fid);
            system('prism mod.sm -mdp -noprobchecks -exporttrans mod.txt');
            fid = fopen('mod.txt','r');
            C1 = textscan(fid,'%f %*f %f %*f %s');
            States=C1{1,1}(1);
            fclose(fid);
            I=[];                          %If deadlock states exist in the flat model,
            for k=1:length(C1{1,1})        %they are shown by self-loops in the flat model
                if isempty(C1{1,3}{k,1})   %with rate 1 with no action name.
                    I=[I,k];               %Since they have no action name, we need to remove them.
                end
            end
            C1{1,1}(I)=[];C1{1,2}(I)=[];C1{1,3}(I)=[];
            C1{1,1}=C1{1,1}+1;C1{1,2}=C1{1,2}+1;
            %%%Tuple
            system('prism mod.sm -exporttransdotstates mod.txt');
            fid = fopen('mod.txt', 'r');      % Open source file.
            prism_tuple = fileread("mod.txt");
            number_of_lines=length(regexp(prism_tuple,'\n','match'));%Total number of lines in the prism_tuple.txt file.
            for k=1:number_of_lines-States-1 % Remove first few lines.
                fgetl(fid) ;
            end
            buffer = fread(fid, Inf) ;                % Read rest of the file.
            fclose(fid);
            fid = fopen('prism_tuple2.txt', 'w');     % Open destination file.
            fwrite(fid, buffer) ;                     % Save to file.
            fclose(fid) ;
            fid = fopen('prism_tuple2.txt', 'r');     % Extract tuples
            clear Tuple2
            for k=1:States
                line=fgetl(fid);
                [~,f1]=regexp(line,'(');
                [~,f2]=regexp(line,')');
                [~,f3]=regexp(line,'[');
                line2=line(f1+1:f2-1);
                line3=line(1:f3-1);
                Tuple2(k,:)=[cell2mat(textscan(line3,'%f')) cell2mat(textscan(line2,...
                    '%f','Delimiter',','))'];
            end
            Tuple2(:,1)=Tuple2(:,1)+1;% In prism, state numbering starts from 0, but I'd like to start from 1.
            clear Anew
            for k=1:length(C1{1,1})
                Anew.matrix(k,:)=[Tuple2(C1{1,1}(k),2:end) , Tuple2(C1{1,2}(k),2:end)];
                Anew.action(k,1)=C1{1,3}(k);
            end
            %%%
            for k=1:length(Tuple)
                first=Tuple(k,im_children(leaf)+1);
                B3=any(Rate.(char(modules(im_children(leaf)))).(char(A.action(ref_t)))(first,:));
                second=Tuple(k,children+1);
                ind=find(ismember(Anew.matrix(:,1:length(children)),second,'rows'));
                B4=any(ismember(Anew.action(ind),A.action(ref_t)));
                if B3 && B4
                    flag=false;
                    break
                end
            end
        else % X1 and X2 are not leaves.
            %X1, X2
            keySet = [1 ,2];
            valueSet1 = nonzeros(cell2mat(alls.child(im_children(1),[3:end])));
            valueSet1 = valueSet1(valueSet1<=number_of_modules);
            valueSet2 = nonzeros(cell2mat(alls.child(im_children(2),[3:end])));
            valueSet2 = valueSet2(valueSet2<=number_of_modules);
            valueSet ={valueSet1,valueSet2};
            M = containers.Map(keySet,valueSet,'UniformValues',false);
            for h=1:2 %left and right children
                children=M(h);
                comb=part1 + newline;
                for k=1:length(children)
                    comb=comb + newline + Single_modul_PRIME(children(k));
                end
                comb=comb + newline + 'system' + newline + syncs{im_children(h)-number_of_modules,7} + newline + 'endsystem';
                fid = fopen('mod.sm','wt');
                fprintf(fid, comb);
                fclose(fid);
                system('prism mod.sm -mdp -noprobchecks -exporttrans mod.txt');
                fid = fopen('mod.txt','r');
                C1 = textscan(fid,'%f %*f %f %*f %s');
                States=C1{1,1}(1);
                fclose(fid);
                I=[];                          %If deadlock states exist in the flat model,
                for k=1:length(C1{1,1})        %they are shown by self-loops in the flat model
                    if isempty(C1{1,3}{k,1})   %with rate 1 with no action name.
                        I=[I,k];               %Since they have no action name, we need to remove them.
                    end
                end
                C1{1,1}(I)=[];C1{1,2}(I)=[];C1{1,3}(I)=[];
                C1{1,1}=C1{1,1}+1;C1{1,2}=C1{1,2}+1;
                %%%Tuple
                system('prism mod.sm -exporttransdotstates mod.txt');
                fid = fopen('mod.txt', 'r');      % Open source file.
                prism_tuple = fileread("mod.txt");
                %length(regexp(prism_tuple,'\\n','match'))
                number_of_lines=length(regexp(prism_tuple,'\n','match'));%Total number of lines in the prism_tuple.txt file.
                for k=1:number_of_lines-States-1 % Remove first few lines.
                    fgetl(fid) ;
                end
                buffer = fread(fid, Inf) ;                % Read rest of the file.
                fclose(fid);
                fid = fopen('prism_tuple2.txt', 'w');     % Open destination file.
                fwrite(fid, buffer) ;                     % Save to file.
                fclose(fid) ;
                fid = fopen('prism_tuple2.txt', 'r');     % Extract tuples
                clear Tuple2
                for k=1:States
                    line=fgetl(fid);
                    [~,f1]=regexp(line,'(');
                    [~,f2]=regexp(line,')');
                    [~,f3]=regexp(line,'[');
                    line2=line(f1+1:f2-1);
                    line3=line(1:f3-1);
                    Tuple2(k,:)=[cell2mat(textscan(line3,'%f')) cell2mat(textscan(line2,...
                        '%f','Delimiter',','))'];
                end
                Tuple2(:,1)=Tuple2(:,1)+1;% In prism, state numbering starts from 0, but I'd like to start from 1.
                clear Anew
                for k=1:length(C1{1,1})
                    Anew.matrix(k,:)=[Tuple2(C1{1,1}(k),2:end) , Tuple2(C1{1,2}(k),2:end)];
                    Anew.action(k,1)=C1{1,3}(k);
                end
                valueSet_Anew{h} =Anew;
            end
            M_Anew = containers.Map(keySet,valueSet_Anew,'UniformValues',false);
            for k=1:length(Tuple)
                first=Tuple(k,M(1)+1);
                ind1=find(ismember(M_Anew(1).matrix(:,1:length(M(1))),first,'rows'));
                B3=any(ismember(M_Anew(1).action(ind1),A.action(ref_t)));
                second=Tuple(k,M(2)+1);
                ind2=find(ismember(M_Anew(2).matrix(:,1:length(M(2))),second,'rows'));
                B4=any(ismember(M_Anew(2).action(ind2),A.action(ref_t)));
                if B3 && B4
                    flag=false;
                    break
                end
            end
%             if ~(B3 && B4)
%                 add_action=[add_action, var2(j)];
%             end
            %clear first second ind1 ind2 M valueSet1 valueSet2 valueSet C1 Tuple2 Buffer States M_Anew B1 B2 B3 B4
        end
        if flag % no type A sp. tr., then check for type B.
            new_syncs=syncs;
            temp=char(new_syncs(var2(j)-number_of_modules,1));
            temp=strcat('|[',A.action(ref_t),',',temp(3:end));
            new_syncs(var2(j)-number_of_modules,1)=temp;
            new_mo=cell(1,1);
            for h2=1:2*size(new_syncs,1)+1
                if ~mod(h2,2)%even
                    new_mo=strcat(new_mo,new_syncs(h2/2));
                else%odd
                    new_mo=strcat(new_mo,SYS_Parts(h2/2+1/2));
                end
            end
            system_modified=new_mo;
            %
            New_Rate=Rate;
            New_actions_struc=actions_struc;
            New_actions_struc.(char(A.action(ref_t)))=[New_actions_struc.(char(A.action(ref_t))), var2(j)];
            COMB=SELF_ALGO(var2(j),alls,number_of_modules,A,ref_t,New_actions_struc);
            Px=nonzeros(cell2mat(alls.child(var2(j),3:end)));
            Px=sort(Px(Px<=number_of_modules));%seq. components under node X->var2(j)
            %%
%             if length(add_action)
%             New_act=act;
%             New_alls=alls;
%             Script_New_must_cannot_may;
%             Find_New_PS;
            %%
            Tx=[];
            for k=1:length(indice2)
                if ~isempty(intersect(Px,cell2mat(A.PS(indice2(k)))))
                    Tx=[Tx,indice2(k)];
                end
            end
            for k=1:size(COMB,1)%for each combination...
                COMBk=nonzeros(COMB(k,:));
                add_selfloop=[];
                for h1=1:length(Tx)
                    SLcandid=setdiff(COMBk,cell2mat(A.PS(Tx(h1))));
                    for h2=1:length(SLcandid)
                        state_with_added_selfloops=A.matrix(Tx(h1),SLcandid(h2));
                        if Rate.(char(modules(SLcandid(h2)))).(A.action{ref_t})...
                                (state_with_added_selfloops,state_with_added_selfloops)
                            continue %If there already exists a selfloop in this state, skip the state and iteration.
                        end
                        if ~isempty(add_selfloop)
                            if ismember([SLcandid(h2) state_with_added_selfloops],add_selfloop,'rows')
                                continue
                            end
                        end
                        add_selfloop=[add_selfloop ; [SLcandid(h2) state_with_added_selfloops]];
                    end
                end
                %add_selfloop=[2 1;2 2];
                if isempty(add_selfloop)
                    continue
                end
                Mod_Single_modul_PRIME=Single_modul_PRIME;
                for h2=1:size(add_selfloop,1)
                    [~,number1]=regexp(Single_modul_PRIME(add_selfloop(h2,1)),'\n','match');
                    state_name=Single_modul_PRIME(add_selfloop(h2,1));
                    state_name=char(state_name);
                    state_name=state_name(number1(1)+1:number1(2)-1);
                    [~,number2]=regexp(state_name,':','match');
                    state_name=state_name(1:number2(1)-1);
                    state_name= state_name(find(~isspace(state_name)));
                    module_with_added_selfloops=char(Mod_Single_modul_PRIME(add_selfloop(h2,1)));
                    self=join(['[',A.action(ref_t),'] ',state_name,'=',string(add_selfloop(h2,2))...
                        ,' -> 1 : (',state_name,'''=',string(add_selfloop(h2,2)),');'],'');
                    New_Rate.(char(modules(add_selfloop(h2,1)))).(char(A.action(ref_t)))(add_selfloop(h2,2),add_selfloop(h2,2))=1;
                    self=strjoin (self,'\n');
                    module_with_added_selfloops=strjoin([module_with_added_selfloops(1:number1(2)+1),self,module_with_added_selfloops(number1(2)+1:end)],'\n');
                    Mod_Single_modul_PRIME(add_selfloop(h2,1))=module_with_added_selfloops;
                end
                new_model=join(Mod_Single_modul_PRIME,newline);
                new_model=strcat(part1,new_model);
                new_model=new_model + newline + 'system' + newline + system_modified + newline + 'endsystem';
                fid = fopen('mod.sm','wt');
                fprintf(fid, new_model);
                fclose(fid);
                system('prism mod.sm -mdp -noprobchecks -exporttrans mod.txt');
                fid = fopen('mod.txt', 'r');
                out=fgetl(fid);
                C = textscan(fid,'%f %*f %f %f %s');
                fclose(fid);
                new=cell2mat(textscan(out,'%n %*n %*n'));
                a=C{1,1}+1;
                b=C{1,2}+1;
                %Trans_New=[a,b];

                for h2=1:length(a)
                    Trans_New(h2,:)={num2str(a(h2)) num2str(b(h2)) char(C{1,4}(h2))} ;
                end
                % Trans_New=cellfun(@num2str,Trans_New,'un',0);
                Trans_New=categorical(Trans_New);
                
                Tran_New_Save{1,1}=[a,b,C{1,3}];% Save before unique command, so the parallel transitions are still in Tran_New_Save
                Tran_New_Save{1,2}=C{1,4};
                I=[];
                for h2=1:length(a)
                    if a(h2)==b(h2)
                        I=[I,h2];
                    end
                end
                   
                
                Trans_New(I,:)=[];% Tran_New shows the unique non_selfloop trs.
                Tran_New_Save{1,1}(I,:)=[];
                Tran_New_Save{1,2}(I,:)=[];
                %Tran_New_Save=Trans_New;
                Trans_New=unique(Trans_New,'rows','stable');
                new(2)=length(Trans_New);
                old=[FlatModelState,FlatNonSelfloopUniqueTrs];
%                 Trans_New=sortrows(Trans_New);
%                 Trans_Old=sortrows(Trans_Old);
                if isequal(new,old) && isequal(Trans_New,Trans_Old)%i.e. successful combination, no type A nor B.
                    %update modules
                    Single_modul_PRIME=Mod_Single_modul_PRIME;
                    %update system
                    temp=char(syncs(var2(j)-number_of_modules,1));
                    temp=strcat('|[',A.action(ref_t),',',temp(3:end));
                    syncs(var2(j)-number_of_modules,1)=temp;
                    new_mo=cell(1,1);
                    for h2=1:2*size(syncs,1)+1
                        if ~mod(h2,2)%even
                            new_mo=strcat(new_mo,syncs(h2/2));
                        else%odd
                            new_mo=strcat(new_mo,SYS_Parts(h2/2+1/2));
                        end
                    end
                    main_system=new_mo;
                    %update Rate matrix and add selfloops.
                    for h2=1:size(add_selfloop,1)
                        Rate.(char(modules(add_selfloop(h2,1)))).(char(A.action(ref_t)))(add_selfloop(h2,2),add_selfloop(h2,2))=1;
                    end
                    %update action_struct
                    actions_struc.(char(A.action(ref_t)))=[actions_struc.(char(A.action(ref_t))) , var2(j)];
                    %update add_action(used in Script_New_must_cannot_may.m)
                    add_action=[add_action, var2(j)];
                    break%-> a successful combination is found, so no need to check the rest of the combinations.
                elseif k==size(COMB,1)
                    %all combinations are checked and no successful one, so
                    %break to the next node X.
                    flag2=false;
                end
            end    
        end
    end%end of 'var2' for loop
    %% Script to find the new must may cannot neighborhoods
    Script_New_must_cannot_may;
    %% Script to check all transitions and find PS for the new system, create the new set of equations.
    New_IS_ref=IS_ref;
    Script_CheckTransitionsNewPS_NewSetOfEq;
    %New_Eq_string=New_Eq_string';
    clc;

    RIS_ref=[];
    for j=1:length(New_IS_ref)
        if ~any(Rate.(char(modules(New_IS_ref(j)))).(char(A.action(ref_t)))(:))
            RIS_ref=[RIS_ref, New_IS_ref(j)];
        end
    end

    
    %%%
    nochange=setdiff(1:number_of_modules,RIS_ref);
    new_model=strjoin(Single_modul(nochange),'\n');
    new_model=strcat(part1,new_model);
    New_Rate=Rate;
    New_actions_struc=actions_struc;
    New_actions_struc.(char(A.action(ref_t)))=[New_actions_struc.(char(A.action(ref_t))), add_action];
    for j=1:length(RIS_ref)
        states_with_added_selfloops=unique(A.matrix(indice2,RIS_ref(j)));
        [~,number1]=regexp(Single_modul(RIS_ref(j)),'\n','match');
        state_name=Single_modul(RIS_ref(j));
        state_name=char(state_name);
        state_name=state_name(number1(1)+1:number1(2)-1);
        [~,number2]=regexp(state_name,':','match');
        state_name=state_name(1:number2(1)-1);
        state_name= state_name(find(~isspace(state_name)));
        module_with_added_selfloops=char(Single_modul(RIS_ref(j)));
        for k=1:length(states_with_added_selfloops)
            self(k)=join(['[',A.action(ref_t),'] ',state_name,'=',string(states_with_added_selfloops(k))...
                ,' -> 1 : (',state_name,'''=',string(states_with_added_selfloops(k)),');'],'');
            New_Rate.(char(modules(RIS_ref(j)))).(char(A.action(ref_t)))(states_with_added_selfloops(k),states_with_added_selfloops(k))=1;
        end
        self=strjoin (self,'\n');
        module_with_added_selfloops=strjoin([module_with_added_selfloops(1:number1(2)+1),self,module_with_added_selfloops(number1(2)+1:end)],'\n');
        new_model={char(new_model),char(module_with_added_selfloops)};
        new_model=strjoin(new_model,'\n');
    end
    new_model=strjoin([new_model,' system',system_modified,'endsystem'],'\n');%This is the new model with added selfloops.
    %% We find 'NewFlatModelTran' and 'NewFlatModelState' and compare
    % them with their original values. If they are not equal, then
    % we are in trouble!
    fid = fopen('New_model.sm','wt');
    fprintf(fid, new_model);
    fclose(fid);
    system('prism New_model.sm -mdp -noprobchecks -exporttrans mod.txt');
    fileID1 = fopen('mod.txt','r');
    C1 = textscan(fileID1,'%f %f %f %f %s');
    len=length(C1{1,1});I=[];     %If deadlock states exist in the flat model,
    for j=2:len                  %they are shown by self-loops in the flat model
        if isempty(C1{1,5}{j,1})   %with rate 1, and they have no action name.
            I=[I,j];             %Since they have no action name, we need to remove them.
        end
    end
    C1{1,1}(I)=[];C1{1,2}(I)=[];
    C1{1,3}(I)=[];C1{1,4}(I)=[];C1{1,5}(I)=[];
    %%%%
    fclose(fileID1);
    NewFlatModelState=C1{1,1}(1);%Number of states in the flat model.
    NewFlatModelTran=length(C1{1,1})-1;
    if NewFlatModelState~=FlatModelState || NewFlatModelTran~=FlatModelTran
        error('After adding selfloops, number of combined model transitions/states is changed, i.e. There are SPURIOUS Transitions.')
    end
    %% Script to find the new must may cannot neighborhoods
    New_IS_ref=IS_ref;
    Script_New_must_cannot_may;
    %% Script to check all transitions and find PS for the new system, create the new set of equations.
    Script_CheckTransitionsNewPS_NewSetOfEq;
    New_Eq_string_in_scope=New_Eq_string;
    solution1_num
    if ~isempty(solution1_num) % A solution is found!
        A.f(~find(rem))=1;
        Tmod=A.f(A.f~=1); %update Tmod
        disp('A solution is found in the original scope. The new system with added selfloops and modified sync sets are saved in "new_model"');
        continue %continue in the main While loop.
    end
   % new_length_add_action=length(add_action);
end
%%
%         add_action=[];
system_modified=char(system_modified);
r=regexp(system_modified,'|[')-1;
l=regexp(system_modified,']|')+1;
r1=regexp(system_modified,'|[')-1;
l1=regexp(system_modified,']|')+1;
for j=1:length(r)
    syncs{j,1}=system_modified(r(j):l(j));
end
for j=1:length(r1)
    syncs{j,5}=[r1(j)];
    syncs{j,6}=[r1(j)+length(syncs{j,1})-1];
end
%%%
%% Arrang paranthesis in matrix "c_mat" in the string "main_system".
[~,right_par_main]=regexp(system_modified,'(','match');
[~,left_par_main]=regexp(system_modified,')','match');
quaran_right=zeros(1,length(right_par_main));
quaran_left=zeros(1,length(left_par_main));
count=1;
while ~all(quaran_right)
    for l=1:length(right_par_main)
        a_temp=right_par_main;
        a_temp(1:l)=[];
        a_temp=a_temp(find(~quaran_right(length(quaran_right)-length(a_temp)+1:end)));
        b=left_par_main(left_par_main>right_par_main(l));
        b=b(find(~quaran_left(length(quaran_left)-length(b)+1:end)));
        if  ~isempty(a_temp) &&  b(1)<a_temp(1) && ~quaran_right(l) && ~quaran_left(find(left_par_main==b(1)))
            %syncs{j,5}=main_system(right_par_main(j):left_par_main(length(quaran_left)-length(b)+1));
            quaran_right(l)=1;
            quaran_left(find(left_par_main==b(1)))=1;
            c_mat(count,[1 2])=[right_par_main(l), b(1)];
            count=count+1;break
        elseif isempty(a_temp)
            quaran_right(l)=1;
            quaran_left(find(left_par_main==b(1)))=1;
            c_mat(count,[1 2])=[right_par_main(l), b(1)];
            count=count+1;break
        end
    end
end
%Create the processes in every sync set.
for l=1:number_of_nodes-number_of_modules
    count=1;temp=[];
    for k=1:size(c_mat,1)
        if  cell2mat(syncs(l,5))>c_mat(k,1)&& cell2mat(syncs(l,5))<c_mat(k,2) && ...
                cell2mat(syncs(l,6))>c_mat(k,1)&& cell2mat(syncs(l,6))<c_mat(k,2)
            temp(count,[1 2])=[c_mat(k,1),c_mat(k,2)];
            count=count+1;
        elseif l==root.indice-number_of_modules
            syncs{l,7}=system_modified;
            count=count+1;
        end
    end
    if ~isempty(temp)
        temp=temp( find(min(abs(temp(:,1)-temp(:,2)))),:);
        syncs{l,7}=system_modified(temp(1):temp(2));
    end
end
%clear children comb a_temp
%%%%------------>var7<---------------
for j=1:length(var7) %Outside the current scope.
    % find the new IS_ref after increasing the scope
    New_IS_ref=cell2mat(alls.child(var7(j),[3:end]));
    New_IS_ref=New_IS_ref(find(New_IS_ref));
    New_IS_ref=New_IS_ref(New_IS_ref<=number_of_modules);
    %%%%%
    children=cell2mat(alls.child(var7(j),3:end));
    children=children(find(children));
    children=children(children<=number_of_modules);
    im_children=cell2mat(alls.child(var7(j),3:4));
    B1=im_children(1)>number_of_modules;
    B2=im_children(2)>number_of_modules;
    new_sync=[syncs{var7(j)-number_of_modules,1}(1:2),char(A.action(ref_t)),',',syncs{var7(j)-number_of_modules,1}(3:end)];
    if B1 && B2
        system_modified=[syncs{im_children(1)-number_of_modules,7},new_sync,syncs{im_children(2)-number_of_modules,7}];
    elseif B1 && ~B2
        system_modified=[syncs{im_children(1)-number_of_modules,7},new_sync,char(modules(im_children(2)))];
    elseif ~B1 && B2
        system_modified=[char(modules(im_children(1))),new_sync,syncs{im_children(2)-number_of_modules,7}];
    else
        system_modified=[char(modules(im_children(1))),new_sync,char(modules(im_children(2)))];
    end
    comb=part1 + newline;
    for k=1:length(children)
        comb=comb + newline + Single_modul(children(k));
    end
    comb=comb + newline + 'system' + newline + system_modified + newline + 'endsystem';
    fid = fopen('mod.sm','wt');
    fprintf(fid, comb);
    fclose(fid);
    system('prism mod.sm -mdp -noprobchecks -exporttrans mod.txt');
    fid = fopen('mod.txt','r');
    C1 = textscan(fid,'%*f %*f %*f %*f %s');
    fclose(fid);
    I=[];                          %If deadlock states exist in the flat model,
    for l=1:length(C1{1,1})        %they are shown by self-loops in the flat model
        if isempty(C1{1,1}{l,1})   %with rate 1 with no action name.
            I=[I,l];               %Since they have no action name, we need to remove them.
        end
    end
    C1{1,1}(I)=[];
    if any(strcmp(A.action{ref_t},C1{1,1}))
        continue
    else
        add_action=[add_action,var7(j)];
    end
    
    %% Modify the "system".
    im=cell2mat(alls.child(var7(j),[3,4]));
    ch_im1=cell2mat(alls.child(im(1),[3:end]));
    ch_im1=[ch_im1 im(1)];
    if isempty(intersect(ch_im1,IS_ref))
        ch_im=im(1);
    end
    ch_im2=cell2mat(alls.child(im(2),[3:end]));
    ch_im2=[ch_im2 im(2)];
    if isempty(intersect(ch_im2,IS_ref))
        ch_im=im(2);
    end
    var2=cell2mat(alls.child(ch_im,[3:end]));
    var2=var2(find(var2)); %#ok<*FNDSB>
    %var2=[Boul_ch,var2];
    var2=var2(var2>number_of_modules);
    var2=intersect(var2,NoActSets);
    for k=1:length(var2)
        %first add the action to one of the sync sets inside the scope.
        current_node=var2(k)-number_of_modules;
        im_children=cell2mat(alls.child(var2(k),[3 4]));%Two immediate children of 'var2(j)'
        B1=im_children(1)<=number_of_modules;%True if it is a leaf
        B2=im_children(2)<=number_of_modules;%True if it is a leaf
        if B1 && B2%Both immediate children are leaves.
            X1=any(Rate.(char(modules(im_children(1)))).(char(A.action(ref_t)))(:));%X1 is true if there is any ref_action exists in the immediate children.
            X2=any(Rate.(char(modules(im_children(2)))).(char(A.action(ref_t)))(:));%X2 is true if there is any ref_action exists in the other immediate children.
            if xor(X1,X2) %if the action exists in one of the immediate children.
                add_action=[add_action, var2(k)];
            end
        elseif xor(B1,B2)%One of the immediate children is leaf.
            if B1
                leaf=1;
                not_leaf=2;
            else
                leaf=2;
                not_leaf=1;
            end
            X1=any(Rate.(char(modules(im_children(leaf)))).(char(A.action(ref_t)))(:));%leaf
            %X2:non_leaf; It is a process.
            children=cell2mat(alls.child(im_children(not_leaf),[3:end]));
            children=children(find(children));
            children=children(children<=number_of_modules);
            comb=part1 + newline;
            for h=1:length(children)
                comb=comb + newline + Single_modul(children(h));
            end
            comb=comb + newline + 'system' + newline + syncs{im_children(not_leaf)-number_of_modules,7} + newline + 'endsystem';
            fid = fopen('mod.sm','wt');
            fprintf(fid, comb);
            fclose(fid);
            system('prism mod.sm -mdp -noprobchecks -exporttrans mod.txt');
            fid = fopen('mod.txt','r');
            C1 = textscan(fid,'%*f %*f %*f %*f %s');
            fclose(fid);
            I=[];                          %If deadlock states exist in the flat model,
            for h=1:length(C1{1,1})        %they are shown by self-loops in the flat model
                if isempty(C1{1,1}{h,1})   %with rate 1 with no action name.
                    I=[I,h];               %Since they have no action name, we need to remove them.
                end
            end
            C1{1,1}(I)=[];X2=false;
            if any(strcmp(A.action{ref_t},C1{1,1}))
                X2=true; %X2 is true so ref_action exists in this process.
            end
            if xor(X1,X2)
                add_action=[add_action, var2(k)];
            end
        else % None of the immediate children is leaf.
            %X1, X2
            for h=1:2
                children=cell2mat(alls.child(im_children(h),[3:end]));
                children=children(find(children));
                children=children(children<=number_of_modules);
                comb=part1 + newline;
                for l=1:length(children)
                    comb=comb + newline + Single_modul(children(l));
                end
                comb=comb + newline + 'system' + newline + syncs{im_children(h)-number_of_modules,7} + newline + 'endsystem';
                fid = fopen('mod.sm','wt');
                fprintf(fid, comb);
                fclose(fid);
                system('prism mod.sm -mdp -noprobchecks -exporttrans mod.txt');
                fid = fopen('mod.txt','r');
                C1 = textscan(fid,'%*f %*f %*f %*f %s');
                fclose(fid);
                I=[];                          %If deadlock states exist in the flat model,
                for l=1:length(C1{1,1})        %they are shown by self-loops in the flat model
                    if isempty(C1{1,1}{l,1})   %with rate 1 with no action name.
                        I=[I,l];               %Since they have no action name, we need to remove them.
                    end
                end
                C1{1,1}(I)=[];X(k)=false;
                if any(strcmp(A.action{ref_t},C1{1,1}))
                    X(k)=true;
                end
            end
            if xor(X(1),X(2))
                add_action=[add_action, var2(k)];
            end
        end
    end%end of 'var2' for loop
    
    %         if isempty(add_action)
    %             disp(['No solution is found for action ',A.action{ref_t}])
    %             continue %continue with the next iteration in the while loop.
    %         end
    %% RIS_ref: Restricted IS_ref
    %RIS_ref=[];
    children=cell2mat(alls.child(ch_im,3:end));
    children=[children ch_im];
    children=children(children<=number_of_modules);
    for k=1:length(children)
        if ~any(Rate.(char(modules(children(k)))).(char(A.action(ref_t)))(:))
            RIS_ref=[RIS_ref, New_IS_ref(k)];
        end
    end
    if isempty(RIS_ref)
        disp(['No solution is found for action ',A.action{ref_t}])
        continue %continue with the next iteration in the var7 for loop.
    end
    
    
    %% Add self-loops. Also in this part the "NewRate" is created which is the transition rate matrices of all the components after adding the selfloops.
    nochange=setdiff(1:number_of_modules,RIS_ref);
    new_model=strjoin(Single_modul(nochange),'\n');
    new_model=strcat(part1,new_model);
    New_Rate=Rate;
    New_actions_struc=actions_struc;
    New_actions_struc.(char(A.action(ref_t)))=[New_actions_struc.(char(A.action(ref_t))), add_action];
    for k=1:length(RIS_ref)
        states_with_added_selfloops=unique(A.matrix(indice2,RIS_ref(k)));
        [~,number1]=regexp(Single_modul(RIS_ref(k)),'\n','match');
        state_name=Single_modul(RIS_ref(k));
        state_name=char(state_name);
        state_name=state_name(number1(1)+1:number1(2)-1);
        [~,number2]=regexp(state_name,':','match');
        state_name=state_name(1:number2(1)-1);
        state_name= state_name(find(~isspace(state_name)));
        module_with_added_selfloops=char(Single_modul(RIS_ref(k)));
        for h=1:length(states_with_added_selfloops)
            self(h)=join(['[',A.action(ref_t),'] ',state_name,'=',string(states_with_added_selfloops(h))...
                ,' -> 1 : (',state_name,'''=',string(states_with_added_selfloops(h)),');'],'');
            New_Rate.(char(modules(RIS_ref(k)))).(char(A.action(ref_t)))(states_with_added_selfloops(h),states_with_added_selfloops(h))=1;
        end
        self=strjoin (self,'\n');
        module_with_added_selfloops=strjoin([module_with_added_selfloops(1:number1(2)+1),self,module_with_added_selfloops(number1(2)+1:end)],'\n');
        new_model={char(new_model),char(module_with_added_selfloops)};
        new_model=strjoin(new_model,'\n');
    end
    system_modified=string(system_modified);
    new_model=strjoin([new_model,' system',system_modified,'endsystem'],'\n');%This is the new model with added selfloops.
    %% We find 'NewFlatModelTran' and 'NewFlatModelState' and compare
    % them with their original values. If they are not equal, then
    % we are in trouble!
    fid = fopen('mod.sm','wt');
    fprintf(fid, new_model);
    fclose(fid);
    system('prism mod.sm -mdp -noprobchecks -exporttrans mod.txt');
    fileID1 = fopen('mod.txt','r');
    C1 = textscan(fileID1,'%f %f %f %f %s');
    len=length(C1{1,1});I=[];     %If deadlock states exist in the flat model,
    for k=2:len                  %they are shown by self-loops in the flat model
        if isempty(C1{1,5}{k,1})   %with rate 1, and they have no action name.
            I=[I,k];             %Since they have no action name, we need to remove them.
        end
    end
    C1{1,1}(I)=[];C1{1,2}(I)=[];
    C1{1,3}(I)=[];C1{1,4}(I)=[];C1{1,5}(I)=[];
    %%%%
    fclose(fileID1);
    NewFlatModelState=C1{1,1}(1);%Number of states in the flat model.
    NewFlatModelTran=length(C1{1,1})-1;
    if NewFlatModelState~=FlatModelState || NewFlatModelTran~=FlatModelTran
        error('After adding selfloops, number of combined model transitions/states is changed, i.e. There are SPURIOUS Transitions.')
    end
    %% Script to find the new must may cannot neighborhoods
    Script_New_must_cannot_may;
    %% Script to check all transitions and find PS for the new system, create the new set of equations.
    Script_CheckTransitionsNewPS_NewSetOfEq;
    New_Eq_string=New_Eq_string';
    if ~isempty(solution1_num) % A solution is found!
        A.f(find(~rem))=1;
        Tmod=A.f(A.f~=1); %update Tmod
        disp('A solution is found in the original scope. The new system with added selfloops and modified sync sets are saved in "new_model"');
        break %continue in the main While loop.
    end
end%end of var7 for loop.
%%%%%%%


%end %main while loop
% end




