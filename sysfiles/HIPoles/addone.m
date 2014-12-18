load StartingPoles.dat

StartingPoles=[StartingPoles(1) ; StartingPoles];
fid = fopen('StartingPoles.dat','w');
fprintf(fid,'%E\n',StartingPoles);
fclose(fid);
