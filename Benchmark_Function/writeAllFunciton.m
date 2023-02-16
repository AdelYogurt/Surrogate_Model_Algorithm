clc;
clear;
close all hidden;

write_file=fopen('BenchmarkFunction.m','w');

all_files=dir(pwd);

for file_index=1:length(all_files)
    file_name=all_files(file_index).name;
    if length(file_name) > 8 && strcmp(file_name(1:8),'function')
        read_file=fopen(file_name,'r');
        while ~feof(read_file)
            string=fgetl(read_file);
            fprintf(write_file,[string,'\n']);
        end
        fclose(read_file);
        clear("read_file");
    end
end


fclose(write_file);
clear("write_file");
