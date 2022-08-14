%%Right Hand Side of the equations
Pos=cell2mat(regexp(New_Eq_string,'=='));
for i=1: length(New_Eq_string)
    temp=char(New_Eq_string(i));
    rhs(i)=eval(temp(Pos(i)+2:end));
end
rhs_NEW=rhs;
clear temp
variables_all=[];
for i=1: length(New_Eq_string)
    temp=symvar(New_Eq_string(i));
    variables_all=[variables_all;temp];
    variables_all=unique(variables_all);
end
clear temp
keySet = variables_all;
ValueSet_logic=ones(size(variables_all),'logical');
ValueSet_Value=zeros(size(variables_all));
M_logic = containers.Map(keySet,ValueSet_logic);
M_Value = containers.Map(keySet,ValueSet_Value);
while any(cell2mat(values(M_logic))) 
    for i=1: length(New_Eq_string)
        eq=New_Eq_string(i);
        variables=symvar(eq);
        if i==1 % The first equation
            for j=1:length(variables)-1
                M_Value(char(variables(j)))=1;
                M_logic(char(variables(j)))=false;
            end
            M_Value(char(variables(length(variables))))=rhs(1);
            M_logic(char(variables(length(variables))))=false;
        else % The rest of the equations
            Prod=1;
            for j=1:length(variables)
                temp(j)=M_logic(char(variables(j)));
                if ~temp(j)
                    Prod=M_Value(char(variables(j)))*Prod;
                end
            end
            temp2=find(temp==true);
            if length(temp2)==1
                M_Value(char(variables(temp2)))=rhs(i)/Prod;
                M_logic(char(variables(temp2)))=false;
            elseif isempty(temp2) 
                if  Prod==rhs(i)
                    continue
                else
                    rhs_NEW(i)=Prod;
                end
            else
                error('Nashod!')
            end
        end 
    end
end
file_name='Satisfiying_f_Poll_'+string(number_of_modules-1)+'.mat';
save( file_name , 'rhs_NEW')

Modified_New_Eq_string=New_Eq_string;
for i=1: length(New_Eq_string)
    temp=char(New_Eq_string(i));
    temp=[temp(1:Pos(i)+1),num2str(rhs_NEW(i))];
    Modified_New_Eq_string(i)=temp;
end
equations=[];
for i=1:length(Modified_New_Eq_string)
    equations=[equations;eval(Modified_New_Eq_string(i))];
end
assume(symvar(equations)>0);% assume that all the variables are positive and non-zero.
solution1=solve(equations,symvar(equations), 'IgnoreAnalyticConstraints',1);

file_name='Satisfiying_f_Poll_'+string(number_of_modules-1)+'.mat';
save( file_name , 'rhs_NEW') 
