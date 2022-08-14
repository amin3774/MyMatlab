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