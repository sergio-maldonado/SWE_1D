function printfile(head,foot,filename,format,Wr,path)
% Function that prints a matrix into a file

X = sprintf('Writing file %s ...',filename);
disp(X)
cd(path)

fid = fopen(filename,'w');
fprintf(fid, head);
fprintf(fid,format,Wr);
fprintf(fid, foot);
fclose(fid);

disp('File written.')