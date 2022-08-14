function COMBS=SELF_ALGO(CH,alls,number_of_modules,A,ref_t,actions_struc)
ch=nonzeros(cell2mat(alls.child(CH,3:end)));
ch=[ch;CH];
remain=ch;
len=length(ch);
comb=cell(1,len);
while ~isempty(remain)
    for j=1:len % check the leafs.
        if ch(j)<=number_of_modules
            comb{j}=ch(j);
            remain(remain==ch(j))=[];
        end
    end
    for j=1:len % check the nodes.
        if ch(j)>number_of_modules && ~ismember(alls.child{ch(j),3},remain) &&...
                ~ismember(alls.child{ch(j),4},remain) && ismember(ch(j),remain)
            left=find(ch==alls.child{ch(j),3});
            right=find(ch==alls.child{ch(j),4});
            if ismember(ch(j),actions_struc.(char(A.action(ref_t))))
                mat1=comb{left};
                mat2=comb{right};
                ma=size(mat1,1);
                mb=size(mat2,1);
                [b1,b2]=ndgrid(1:ma,1:mb);
                comb{j}=[mat1(b1,:),mat2(b2,:)];
                remain(remain==ch(j))=[];
            else %Union of left and right leaves.
                sx = size(comb{left});
                sy = size(comb{right});
                am = max(sx(2),sy(2));
                comb{j} = [[comb{left},zeros(abs([0 am]-sx))];[comb{right},zeros(abs([0,am]-sy))]];
                remain(remain==ch(j))=[];
            end
        end
    end
end
COMBS=cell2mat(comb(ch==CH));