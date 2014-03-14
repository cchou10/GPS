% combo.m		(actual file name: combo.m)
%
%  Copyright (c) 2008 Paul M. Kinter, Jr., and 
%  Copyright (c) 2012 Mark L. Psiaki.  All rights reserved. 
%
% calcultates all possible combinations of length
% 'v_size' using elements in 'vector'
%
% input: 'vector' individual elements of the combinations
% 			'v_size' number of elements in each combination
%					[ vector, v_size ]
%
function combination = combo(vector,v_size)
combination=[];
% if the size of each combination equals the size of the input vector, there
% is only 1 combination
if v_size==1
  combination=[vector]';
else
  % form an n by m matrix where n is the number of combinations and m is
  % the number of elements in each combination
  array=zeros(nchoosek(length(vector),v_size),v_size);
  num_uses=[];
  num_comb=[];
  % fill in the matrix with all combinations
  for i=length(vector):(-1):v_size
    num_comb=[num_comb,nchoosek(i,v_size)];
  end
  for i=1:(length(num_comb)-1)
    num_uses=[num_uses,num_comb(i)-num_comb(i+1)];
  end
  num_uses=[num_uses,1];
  rowA=1;
  rowB=1;
  for i=1:length(num_uses)
    for j=1:num_uses(i)
      array(rowA,1)=vector(i);
      rowA=rowA+1;
    end
    vector2=[];
    for k=i+1:length(vector)
      vector2=[vector2,vector(k)];
    end
    array(rowB:rowB+(num_uses(i)-1),2:size(array,2))=combo(vector2,v_size-1);
    rowB=rowB+num_uses(i);
  end
  combination=array;
end
return
