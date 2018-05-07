%Function to convert Matlab matrix to .dxf
%Input is filename (something.dxf), X vector, Y vector, Z vector
%Input 'clic' is a vector containing number of points for each polyline
%
% %example: clic = [{3},{5}]
% %meaning : the first polyline has 3 points, the second one has 5
%
%make sure that the coordinate vectors have an extra point which is the
%same as the first point, in order to close the polyline.
%
% %example for 3 point polyline: X = [234.123; 456.745; 342.230; 234.123]
% %notice that although the polyline has only 3 points, there are 4 points
% %in the vector. Notice also that the last point is the same as the 1st
%
%based on the function writedxf.m from Greg Siegle. Refer to Matlab 
%File Exchange, search for writedxf.m for more info.
%
%Written by Arnadi MURTIYOSO
%ifp Stuttgart & INSA de Strasbourg
%last update: 25/06/15

function[]=polydxf(fname,X,Y,Z,clic) 

fullname=sprintf('%s',fname);
fid=fopen(fullname,'w');
e=1;
f=1;
clic = cell2mat(clic); %turn off if clic is not a cell array
[r,~] = size(clic);
for b=1:r
    a = clic (b,1);
    fprintf(fid,'0\nSECTION\n2\nENTITIES\n');
    for c = e:(e+a)    
    if c < (e+a)
      fprintf(fid,'0\nLINE\n8\n0\n');
      %create new line element
      fprintf(fid,'10\n%.4f\n20\n%.4f\n30\n%.4f\n',X(c),Y(c),Z(c));
      %first coordinate triple - starting point of line-element
      fprintf(fid,'11\n%.4f\n21\n%.4f\n31\n%.4f\n',X(c+1),Y(c+1),Z(c+1));      
      %second coordinate triple - ending point of line-element
    else
        break
    end
    end
    fprintf(fid,'0\nENDSEC\n');
    e = e+clic(b,1)+1;
    f = f+1;
end
fprintf(fid,'0\nEOF\n');
fclose(fid);

