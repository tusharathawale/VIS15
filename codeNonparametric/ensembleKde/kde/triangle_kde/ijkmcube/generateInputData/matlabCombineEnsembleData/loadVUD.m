function [volume, lattice, datatype, scaling, offset] = loadVUD(fileName)
% todo, add more parameters that are returned, so that a file can
% be loaded and written verbatim without any additional work
%
%
fid = fopen(fileName);
aline = fgetl(fid); % vu data version, a comment
aline = fgetl(fid); % Volume Data
aline = fgetl(fid); % Binary

aline = fgetl(fid); % lattice
lattice_token = sscanf(aline, 'DATASET %s');
lattice=strtok(lattice_token, '_');

aline = fgetl(fid); % unimodal

HeaderLine = fgetl(fid); % x y z dims
dimensions = sscanf(HeaderLine, 'DIMENSIONS %d %d %d %d');
xsize = dimensions(1);
ysize = dimensions(2);
zsize = dimensions(3);

volume = zeros(xsize, ysize, zsize);
aline = fgetl(fid); % origin, which is always 0

aline = fgetl(fid); % offset
offset = sscanf(aline, 'OFFSET %f %f %f %f')';

aline = fgetl(fid); % scaling
scaling = sscanf(aline, 'SCALING %f %f %f %f')';

aline = fgetl(fid); % spacing
aline = fgetl(fid); % POINT_DATA num data points
aline = fgetl(fid); % datatype: byte or float

datatype = sscanf(aline, 'SCALARS data %s');
disp(sprintf('datatype is %s', datatype));

aline = fgetl(fid);
if(~isequal(aline,'LOOKUP_TABLE default'))
  disp(sprintf('error in the headerfile, you need a redundant line at the end of the header, but %s\n', aline));
  return;
end

if(isequal(datatype, 'byte'))
  read_datatype = 'uint8';
elseif(isequal(datatype, 'float'))
  read_datatype = 'float32';
else
  disp('ERROR IN RECOVERING FILE DATATYPE');
end

for i = 1:zsize
  [volume(:,:,i), count] = fread(fid, [xsize ysize], read_datatype);
end

fclose(fid);
