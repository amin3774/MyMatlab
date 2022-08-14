counter=1;eqns=[];rem=ones(length(A.f),1);COUN=1;
%%
%% Create Structure A
A_new.matrix(:,[1:number_of_modules])=Tuple(Tran_New_Save{1,1}(:,1),[2:end]);
A_new.matrix(:,[number_of_modules+1:2*number_of_modules])=Tuple(Tran_New_Save{1,1}(:,2),[2:end]);
A_new.action=Tran_New_Save{1,2};
A_new.rate=Tran_New_Save{1,1}(:,3);
%%
indice2_new=find(strcmp(A_new.action,A.action{ref_t,1}));


    for i1=1:length(indice2)
        if ~rem(indice2(i1))
            continue
        end
        counter2=1;
        MS=find(A.matrix(indice2(i1),1:number_of_modules)~=A.matrix(indice2(i1),...
            number_of_modules+1:end));      %Moving Set of t.
        if isempty(MS) % If there is a selfloop in the combined model, we ignore it.
            continue
        end
        SS=setdiff(1:number_of_modules,MS); %SS is the complement of MS.
        %%
        PS=MS; %Initialise
        PS_candid=[];
        % To find PS.
        for h1=1:length(SS)
            for k1=1:length(MS)
                Pi=SS(h1);Pj=MS(k1);
                if ismember(Pi,New_act.(char(A.action(ref_t))).(New_alls.child{Pj,1}).MUST)
                    PS=[PS Pi];
                end
                if ismember(Pi,New_act.(char(A.action(ref_t))).(New_alls.child{Pj,1}).MAY)...
                        && New_Rate.(char(modules(Pi))).(char(A.action(ref_t)))(A.matrix(indice2(i1),Pi),A.matrix(indice2(i1),Pi))~=0
                    PS_candid=[PS_candid Pi];
                end
            end
        end
        PS=unique(PS);
        PS_candid=unique(PS_candid);
        cannots=[];
        for h1=1:length(MS)
            cannots=[cannots; New_act.(char(A.action(ref_t))).(New_alls.child{MS(h1),1}).CANNOT];
        end
        PS_candid=setdiff(PS_candid,cannots);
        PS_candid=setdiff(PS_candid,PS);
        %There might be a component which is in a may neighborhood of a
        %moving component AND at the same time in a must neghborhood of another moving component.
        %This component is for sure participating and that is why it is removed from PS_candid.
        %% Check PS_candid components
        for h1=1:length(PS_candid)
            % Find the smallest subtree,"The_node" consisting of PS_candid and a moving
            % component
            Bol=true;
            %"Bol" is a boolean varaible. For every component in PS_candid,
            %if 'Bol' is true at end, then it is a participating component.
            The_node=[];The_node_len=[];
            a1=find(strcmp(actions(:,1),A.action(ref_t)));
            for k1=1:number_of_nodes-number_of_modules
                if ismember(number_of_modules+k1,cell2mat(actions(a1,[2:end]))) %Check the nodes which incldes ref_action.
                    am=cell2mat(New_alls.child(k1+number_of_modules,[3:end]));
                    am=am(find(am));
                    if ismember(PS_candid(h1),am) && ~isempty(intersect(MS,am))
                        The_node=[The_node,k1+number_of_modules];
                        The_node_len=[The_node_len,length(find(cell2mat(New_alls.child(k1+number_of_modules,[3:end]))))];
                    end
                end
            end
            [~,The_node_len]=min(The_node_len);
            The_node=The_node(The_node_len);
            [parent,~]=find(cell2mat(New_alls.child(number_of_modules+1:end,[3,4]))==PS_candid(h1));
            parent=parent+number_of_modules;
            while parent ~= The_node
                if ismember(parent,cell2mat(actions(a1,[2:end])))
                    children=cell2mat(New_alls.child(parent,[3:end]));
                    children=children(find(children));
                    children=children(children<=number_of_modules);
                    comb=part1 + newline;
                    for k1=1:length(children)
                        comb=comb + newline + Single_modul(children(k1));
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
                    for k1=1:number_of_lines-FlatModelState2-1 % Remove first few lines.
                        fgetl(fid) ;
                    end
                    buffer = fread(fid, Inf) ;             % Read rest of the file.
                    fclose(fid);
                    fid = fopen('mod_tuple.txt', 'w');     % Open destination file.
                    fwrite(fid, buffer) ;                  % Save to file.
                    fclose(fid) ;
                    fid = fopen('mod_tuple.txt', 'r');     % Extract tuples
                    for k1=1:FlatModelState2
                        line=fgetl(fid);
                        [~,f1]=regexp(line,'(');
                        [~,f2]=regexp(line,')');
                        [~,f3]=regexp(line,'[');
                        line2=line(f1+1:f2-1);
                        line3=line(1:f3-1);
                        Tuple2(k1,:)=[cell2mat(textscan(line3,'%f')) cell2mat(textscan(line2,...
                            '%f','Delimiter',','))'];
                    end
                    Tuple2(:,1)=Tuple2(:,1)+1;% In prism, state numbering starts from 0, but
                    % I'd like to start from 1.
                    ind=find(ismember( Tuple2(:,2:end),A.matrix(indice2(i1),children),'rows'));
                    if ismember([ind ind],C2.transitions,'rows')
                        ind2=find(ismember(C2.transitions,[ind ind],'rows'));
                        if any(strcmp(C2.action(ind2),A.action(ref_t)))
                            Bol=true;
                        else
                            Bol=false;
                            i1,h1
                            break
                        end
                    else
                        Bol=false;
                        i1,h1
                        break
                    end
                end
                [parent,~]=find(cell2mat(New_alls.child(number_of_modules+1:end,[3,4]))==parent);% The new parent is the parent of the old one.
                parent=parent+number_of_modules;
            end
            if Bol
                PS=[PS, PS_candid(h1)];
            end
        end
        
        %%
        % If PS(t) is not a subset of IS(ref_t),then skip this t and take the
        % next transition of 'indice2'.
        if ~all(ismember(PS, New_IS_ref))
            continue
        else
        rem(indice2(i1))=false;    
            if ismember(3,PS) && ~ismember(3,MS) 
                WH(COUN)=i1;
                COUN=COUN+1;
            end
            PS_t.vector(counter,:)=zeros(1,number_of_modules);
            PS_t.vector(counter,[1:length(PS)])=PS;
            PS_t.indice(counter)=indice2(i1);
    
            %Obtain the 'variable', where the first column shows the indice
            %of the variable, the second one is the indice
            %of the component, third one shows the source state and the
            %fourth one shows the target state in that component.
            
            if ismember(PS,MS) % if PS=MS
                eqn=[];
                for h1=1:length(PS)
                    variable(counter2,:)=[PS(h1),A.matrix(indice2(i1),PS(h1)),...
                        A.matrix(indice2(i1),PS(h1)+number_of_modules)];
                    eqn=[eqn, variable(counter2,:)];
                    vari{h1}=['X',num2str(PS(h1)),'_',num2str(A.matrix(indice2(i1),PS(h1))),...
                        '_',num2str(A.matrix(indice2(i1),PS(h1)+number_of_modules))];
                    syms(vari{h1})
                    counter2=counter2+1;
                end
                e=string(join(vari,'*'));
                New_Eq_string(counter)=strcat(e,'==',num2str(A.rate(indice2(i1))),'*',num2str(A.f(indice2(i1))));
                %                 eqns.lhs(counter,1:length(eqn))=eqn;
                %                 eqns.rhs(counter,:)=A.rate(indice2(i))*A.f(indice2(i));
                counter=counter+1;
            else %(RSLC)If there are some self-loops and |PS|>|MS|.
                %% RSLC: find the smallest subtree consists of PS: variable name:"local_root"
                a=cell2mat(New_alls.child(number_of_modules+1:end,3:end));
                b=[];c=[];
                for h1=1:number_of_nodes-number_of_modules
                    if ismember(PS,a(h1,:))
                        b=[b, number_of_modules+h1];
                        c=[c,length(find(a(h1,:)))];
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
                    for h1=1:len % check the leafs.
                        if ch(1,h1)<=number_of_modules
                            musts=[];
                            for k1=1:length(MS)
                                musts=[musts; New_act.(char(A.action(ref_t))).(New_alls.child{MS(k1),1}).MUST];
                                musts=unique(musts);
                            end
                            if ismember(ch(1,h1),PS) && ~ismember(ch(1,h1),MS) && ~ismember(ch(1,h1),musts)
                                comb{h1}=ch(1,h1);
                                remain(remain==ch(1,h1))=[];
                            else
                                %comb{j}=0;
                                remain(remain==ch(1,h1))=[];
                            end
                        end
                    end
                    for h1=1:len % check the nodes.
                        if ch(1,h1)>number_of_modules && ~ismember(New_alls.child{ch(h1),3},remain) &&...
                                ~ismember(New_alls.child{ch(h1),4},remain) && ismember(ch(h1),remain)
                            left=find(ch==New_alls.child{ch(h1),3});
                            right=find(ch==New_alls.child{ch(h1),4});
                            if ismember(ch(1,h1),New_actions_struc.(char(A.action(ref_t))))
                                if isempty(comb{left}) && isempty(comb{right})
                                    %comb{j}=[];
                                    remain(remain==ch(1,h1))=[];
                                elseif ~isempty(comb{left}) && isempty(comb{right})
                                    comb{h1}=comb{left};
                                    remain(remain==ch(1,h1))=[];
                                elseif isempty(comb{left}) && ~isempty(comb{right})
                                    comb{h1}=comb{right};
                                    remain(remain==ch(1,h1))=[];
                                elseif ~isempty(comb{left}) && ~isempty(comb{right})
                                    mat1=comb{left};
                                    mat2=comb{right};
                                    ma=size(mat1,1);
                                    mb=size(mat2,1);
                                    [b1,b2]=ndgrid(1:ma,1:mb);
                                    comb{h1}=[mat1(b1,:),mat2(b2,:)];
                                    remain(remain==ch(1,h1))=[];
                                end
                            else %Union of left and right leaves.
                                sx = size(comb{left});
                                sy = size(comb{right});
                                am = max(sx(2),sy(2));
                                comb{h1} = [[comb{left},zeros(abs([0 am]-sx))];[comb{right},zeros(abs([0,am]-sy))]];
                                remain(remain==ch(1,h1))=[];
                            end
                        end
                    end
                end
                % Here we create the equations if |PS|>|MS|, and sebsequently there exists
                %rslc set. First, the vaiables corresponding MS, then
                %variables corresponding Q and finally rslc.
                rslc=comb{find(ch==local_root)} ; %RSLC result.
                multi_trans_idx=indice2(find(ismember(A.matrix(indice2,:),A.matrix(indice2(i1),:),'rows')));
%                 if ~isempty(rslc) && size(rslc,1)~=length(multi_trans_idx)
%                     error('Error 1 in RSLC function. The number of self-loop combinations is different with the number of multi transitions.')
%                 end
%                 if isempty(rslc) && length(multi_trans_idx)>1
%                     error('Error 2 in RSLC function. RSLC is empty but there are some multi transitions.')
%                 end
                eqn=[];vari=[];
                for h1=1:length(MS) %counter2 is the number of variables.
                    variable(counter2,:)=[MS(h1),A.matrix(indice2(i1),MS(h1)),...
                        A.matrix(indice2(i1),MS(h1)+number_of_modules)];
                    eqn=[eqn, variable(counter2,:)];
                    vari{h1}=['X',num2str(MS(h1)),'_',num2str(A.matrix(indice2(i1),MS(h1))),...
                        '_',num2str(A.matrix(indice2(i1),MS(h1)+number_of_modules))];
                    syms(vari{h1})
                    counter2=counter2+1;
                end
                % Q is the set of components defined in line 48 of the algorithm.
                Q=intersect(musts,setdiff(PS,MS));
                for h1=1:length(Q)
                    variable(counter2,:)=[Q(h1),A.matrix(indice2(i1),Q(h1)),...
                        A.matrix(indice2(i1),Q(h1)+number_of_modules)];
                    eqn=[eqn, variable(counter2,:)];
                    vari{counter2}=['X',num2str(Q(h1)),'_',num2str(A.matrix(indice2(i1),Q(h1))),...
                        '_',num2str(A.matrix(indice2(i1),Q(h1)+number_of_modules))];
                    syms(vari{counter2})
                    counter2=counter2+1;
                end
                %rslc
                if ~isempty(rslc) %if rslc is not empty.
                    part3={};e={}; part2={};
                    for h1=1:size(rslc,1)
                        s1=rslc(h1,:);
                        s1=s1(find(s1)); %remove zero elements
                        s2=length(s1);
                        for k1=1:s2
                            part2{h1,k1}=['X',num2str(rslc(h1,k1)),'_',num2str(A.matrix(indice2(i1),rslc(h1,k1))),...
                                '_',num2str(A.matrix(indice2(i1),rslc(h1,k1)+number_of_modules))];
                            syms(part2{h1,k1})
                            counter2=counter2+1;
                        end
                        part3(h1,1:length([vari part2(h1,:)]))=[vari part2(h1,:)];
                        e(h1)=join(part3(h1,:),'*');
                    end
                    complete=string(join(e,'+'));
                    RHS=0;
                    for h1=1:size(rslc,1) %Right Hand Side of the equation, which is the product of f and rates.
                        RHS=RHS+A.rate(multi_trans_idx(h1))*A.f(multi_trans_idx(h1));
                    end
                    New_Eq_string(counter)=strcat(complete,'==',num2str(RHS));
                else % if rslc is empty.
                    e=string(join(vari,'*'));
                    New_Eq_string(counter)=strcat(e,'==',num2str(A.rate(indice2(i1))),'*',num2str(A.f(indice2(i1))));
                end
                rem(multi_trans_idx)=false;
                counter=counter+1;
            end
        end
    end
    New_Eq_string=unique(New_Eq_string);
    New_Eq_string=New_Eq_string';
   % clearvars  -except -regexp ^X New_Eq_string number_of_modules
    equations=[];
    for h1=1:length(New_Eq_string)
        equations=[equations;eval(New_Eq_string(h1))];
    end
    assume(symvar(equations)>0);% assume that all the variables are positive and non-zero.
    solution1=solve(equations,symvar(equations), 'IgnoreAnalyticConstraints',1) %Solve system of eqs.
    solution1_num=[];
    if isstruct (solution1)
        solution1_num = struct2cell(solution1);
        solution1_num = cat(2,solution1_num{:});
        solution1_num=solution1_num';
        solution1_num
    end
%     if ~isempty(solution1_num)
%        fprintf('There exist a solution for this set of modification factors by adding the action to the sync. set(s)%g and adding selfloops in component(s) %g, state(s):[ %g ].', add_action, RIS_ref, states_with_added_selfloops); 
%     end