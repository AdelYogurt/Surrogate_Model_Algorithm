clc;
clear
close all hidden;

fval_nomlz_list=[0.2540;10.0000;0.0743;0.2714;0.0026;0.0572]

[quantile,normal_index]=findUnusual(fval_nomlz_list)

function [quantile,normal_index]=findUnusual(data)
% base on Box diagram to search unusual value
%
x_number=length(data);
[data,index_list]=sort(data);

Q1=getQuantile(data,0.25);
Q3=getQuantile(data,0.75);
IQR=Q3-Q1;

normal_index=1:x_number;
normal_index(data < (Q1-1.5*IQR))=[];
normal_index(data > (Q3+1.5*IQR))=[];

quantile=[min(data(normal_index));
    getQuantile(data(normal_index),0.25);
    getQuantile(data(normal_index),0.5);
    getQuantile(data(normal_index),0.75);
    max(data(normal_index))];

[~,index_list]=sort(index_list);
normal_index=index_list(normal_index);

    function quantile=getQuantile(data,percent)
        index=x_number*percent;
        if (index == fix(index))
            % index is integer
            quantile=0.5*(data(index)+data(index+1));
        else
            quantile=0.5*data(ceil(index));
        end
    end
end