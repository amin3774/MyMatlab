counter=1;eqns=[];rem=ones(length(A.f),1);COUN=1;
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
    end