function ima = normalization_UWA(M,N,p)


%   normalization_UWA compute the contribution/weights of each patch of 
%   size p x p in an image of size M x N.
%                 
%   INPUT:
%     M,N             : image size 
%     p               : patch size
%
%   OUTPUT:
%    ima              : normalization, counts the number of patches a
%                       pixel belong to. usually it is p x p, but one
%                       needs to deal with border issues
%
%   ima = normalization_UWA(M,N,p)
%   
%   Copyright (C) 2012 NLSPCA project
%   Joseph Salmon, Zachary Harmany, Charles-Alban Deledalle, Rebecca Willett
%
%   See The GNU Public License (GPL)
%
%---------------------------------------------------------------------
%
%   This file is part of NLSPCA.
%
%   NLSPCA is free software: you can redistribute it and/or modify
%   it under the terms of the GNU General Public License as
%   published by the Free Software Foundation, either version 3 of
%   the License, or (at your option) any later version.
%
%   NLSPCA is distributed in the hope that it will be useful,
%   but WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%   GNU General Public License for more details.
%
%   You should have received a copy of the GNU General Public
%   License along with NLSPCA.  If not, see
%   <http://www.gnu.org/licenses/>.
%   Joseph Salmon, Zachary Harmany, Charles-Alban Deledalle, Rebecca Willett
%
%   See The GNU Public License (GPL)
p1=p-1;
M1=M+p1;
N1=N+p1;
ima = zeros(M+p-1, N+p-1);

for i = 1:p
    xrange = mod(((i-1)+(1:p))-1,M1)+1;
    for j = 1:N
    
        yrange = mod(((j-1)+(1:p))-1,N1)+1;                    
        ima(xrange, yrange) =ima(xrange, yrange) + 1;        
    end
end

%down
for i = (M-p+1):M
    xrange = mod(((i-1)+(1:p))-1,M1)+1;
    for j = 1:N
    
        yrange = mod(((j-1)+(1:p))-1,N1)+1;                    
        ima(xrange, yrange) =ima(xrange, yrange) + 1;        
    end
end


%left
 for j = 1:p
     yrange = mod(((j-1)+(1:p))-1,N1)+1;                    
     for i = (p+1):(M-p)
   
         xrange = mod(((i-1)+(1:p))-1,M1)+1;
        
        ima(xrange, yrange) =ima(xrange, yrange) + 1;        
    end
end



%right

for j = (N-p+1):N
    
    yrange = mod(((j-1)+(1:p))-1,N1)+1;

    
    for i = (p+1):(M-p)
    
        xrange = mod(((i-1)+(1:p))-1,M1)+1;
                    
        ima(xrange, yrange) =ima(xrange, yrange) + 1;        
    end
end

ima(p:(M), p:(N)) = p*p;
