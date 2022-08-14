%% CANNOT for the new model with added selfloops.
            New_act=act;
            New_alls=alls;
            for h1=1:length(add_action)
                temp1=char(New_alls.child(add_action(h1),2));
                temp1=[temp1(1:2),char(A.action(ref_t)),',',temp1(3:end)];
                New_alls.child(add_action(h1),2)=cellstr(temp1);
            end   
            
            current_act=char(A.action(ref_t));
            for h1=1:number_of_modules
                current_mod=string(New_alls.child(h1,1));
                New_act.(current_act).(current_mod).CANNOT=[];
                [parent,~]=find(cell2mat(New_alls.child(number_of_modules+1:end,[3 4]))==h1);
                parent=parent+number_of_modules; %Initialisation
                old_parent=h1;                    %Initialisation
                cond=true;
                while cond
                    salar=char(New_alls.child(parent,2));
                    salar=salar(3:end-2);      % Remove |[ and ]|.
                    salar(isspace(salar)) = [];% Remove space from the string.
                    salar=split(salar,",");    % Split the action names.
                    if isempty(find(strcmp(current_act,salar), 1))% if current_act does Not exist in the parent node...
                        temp=cell2mat(New_alls.child(parent,[3 4]));
                        temp=temp(temp~=old_parent);
                        if temp<=number_of_modules
                            New_act.(current_act).(current_mod).CANNOT=...
                                [New_act.(current_act).(current_mod).CANNOT;temp];
                        else
                            temp2=nonzeros(cell2mat(New_alls.child(temp,3:end)));
                            temp2=temp2(temp2<=number_of_modules);
                            New_act.(current_act).(current_mod).CANNOT=...
                                [New_act.(current_act).(current_mod).CANNOT;temp2];
                        end
                    end
                    if parent==root.indice
                        cond=false;
                    else
                        old_parent=parent;
                        [parent,~]=find(cell2mat(New_alls.child(number_of_modules+1:end,[3 4]))==parent);
                        parent=parent+number_of_modules;
                    end
                end
            end
        
        %% MUST for the new model with added selfloops.
            for h1=1:number_of_modules
                current_mod=string(New_alls.child(h1,1));
                New_act.(current_act).(current_mod).MUST=[];
                [parent,~]=find(cell2mat(New_alls.child(number_of_modules+1:end,[3 4]))==h1);
                parent=parent+number_of_modules; %Initialisation
                old_parent=h1;                    %Initialisation
                condition=true;
                while condition
                    salar=char(New_alls.child(parent,2));
                    salar=salar(3:end-2);      % Remove |[ and ]|.
                    salar(isspace(salar)) = [];% Remove space from the string.
                    salar=split(salar,",");    % Split the action names.
                    if ~isempty(find(strcmp(current_act,salar), 1))% if current_act exists in the parent node, ham mirim bala ham paiin...
                        temp=cell2mat(New_alls.child(parent,[3 4]));
                        temp=temp(temp~=old_parent);
                        if temp<=number_of_modules%ok
                            New_act.(current_act).(current_mod).MUST=...
                                [New_act.(current_act).(current_mod).MUST;temp];
                        else%be paiin
                            children=nonzeros(cell2mat(New_alls.child(temp,3:end)))';
                            children=children(children<=number_of_modules);
                            for h2=1:length(children)
                                cond=true;
                                [sub_parent,~]=find(cell2mat(New_alls.child(number_of_modules+1:end,[3 4]))==children(h2));
                                sub_parent=sub_parent+number_of_modules;
                                while cond
                                    salar=char(New_alls.child(sub_parent,2));
                                    salar=salar(3:end-2);      % Remove |[ and ]|.
                                    salar(isspace(salar)) = [];% Remove space from the string.
                                    salar=split(salar,",");    % Split the action names.
                                    if isempty(find(strcmp(current_act,salar), 1))
                                        %if cellfun(@isempty,regexp(current_act,alls.child(sub_parent,2), 'once'))
                                        break
                                    elseif sub_parent==parent
                                        cond=false;
                                        New_act.(current_act).(current_mod).MUST=...
                                            [New_act.(current_act).(current_mod).MUST;children(h2)];
                                    else
                                        [sub_parent,~]=find(cell2mat(New_alls.child(number_of_modules+1:end,[3 4]))==sub_parent);
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
                        [parent,~]=find(cell2mat(New_alls.child(number_of_modules+1:end,[3 4]))==parent);
                        parent=parent+number_of_modules;
                    else  %We move to the upper node.
                        old_parent=parent;
                        if parent==root.indice
                            condition = false;
                        end
                        [parent,~]=find(cell2mat(New_alls.child(number_of_modules+1:end,[3 4]))==parent);
                        parent=parent+number_of_modules;
                        cond1=true;%?
                    end
                end
            end
        %% May for the new model with added selfloops.
        %If a component is niether MUST nor CANNOT, it is MAY.
            for h1=1:number_of_modules
                current_mod=string(New_alls.child(h1,1));
                New_act.(current_act).(current_mod).MAY=[];
                u=union(New_act.(current_act).(current_mod).CANNOT,...
                    New_act.(current_act).(current_mod).MUST);
                New_act.(current_act).(current_mod).MAY=setdiff([1:h1-1,h1+1:number_of_modules],u)';
            end
            clear temp temp1