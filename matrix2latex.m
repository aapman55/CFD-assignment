%% Assignment for AE4134: CFD I
% Created by:
% Zhi-li Liu 4146557
% Jasper van Wensveen 4142179
%
%matrix2latex Creates the latex code in bmatrix format

function [ output ] = matrix2latex( matrix, fileoutput )
% Checks whether the input matrix is sparse. dlmwrite can not work with
% sparse matrices. So it needs to be converted to a full matrix first. This
% might give problems for huge matrices.
if issparse(matrix)
    matrix = full(matrix);
end

% Create a file for the latex output
fid = fopen(fileoutput,'w+');

% open bmatrix environment
fprintf(fid,'%s\n','\begin{bmatrix}');

% Loop over every row
for i=1:size(matrix,1)
    % Indent each line in the latex code for cleanliness
   str = sprintf('\t');
   
   % Extract all components of the current row.
   current = matrix(i,:);
   
   % Loop over all components in the current row
   for j = 1:length(current)
       % For all components not equal to the last entry
       if (j < length(current))
            str = [str, num2str(current(j)),sprintf('\t'),'&',sprintf('\t')];
       % The last entry
       else
            str = [str, num2str(current(j))];
       end
   end
   
   % When the current row is not the last row
   if (i < size(matrix,1))
       str = [str,'\\' ,sprintf('\n')];
   % When the current row is the last row
   else
       str = [str,sprintf('\n')];
   end
   
   % Write the corresponding line break
   fprintf(fid,'%s',str);
end

% Close bmatrix environment
fprintf(fid,'%s\n','\end{bmatrix}');

% close file
fclose(fid);
end

