function [face_normal,face_normal_unit,face_area,face_centroid]=findTriFaceNormalArea(face_nodes)
%% Function to calculate area and normal of Tri Element
% Created by Thor Andreassen
% 1/8/22
%
%
% the function takes in a set of 3 nodes, that correspond to the tri
% element corner nodes, and calculates the surface normal direction as well
% as the area of the element.
%
%
% Permission is hereby granted, free of charge, to any person obtaining a
% copy of this software and associated documentation files (the
% "Software"), to deal in the Software without restriction, including
% without limitation the rights to use, copy, modify, merge, publish,
% distribute, sublicense, and/or sell copies of the Software, and to permit
% persons to whom the Software is furnished to do so, subject to the
% following conditions:
%
% The above copyright notice and this permission notice shall be included
% in all copies or substantial portions of the Software.
%
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
% OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
% MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN
% NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
% DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR
% OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE
% USE OR OTHER DEALINGS IN THE SOFTWARE.
%
        vec1=face_nodes(2,:)-face_nodes(1,:);
        vec1_norm=vec1/norm(vec1);
        vec2=face_nodes(3,:)-face_nodes(1,:);
        vec2_norm=vec2/norm(vec2);
        face_area=abs((0.5)*det([1,1,1;vec1;vec2]));
        face_normal=cross(vec1_norm,vec2_norm);
        face_normal_unit=face_normal/norm(face_normal);
        face_centroid=mean(face_nodes);
end
