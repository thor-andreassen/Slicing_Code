function [quality_indexes,min_angles]=getMeshQuality(faces,nodes)
        %% Function to calculate the quality of a given triangulated Mesh
        % Created by Thor Andreassen
        % 1/10/22
        %
        % This function takes in a triangulated mesh as a series of faces
        % and nodes, and returns the corresponding quality indexes as the
        % angles of the 3 sides of every face, as well as the minimum angle
        % of the entire face.
        %
        % Inputs:
                % faces: (n x 3) matrix of connections of nodes defined by
                        % corner node indices for a triangulated mesh
                % nodes: (mx3) matrix of node coordinates (x, y, z) with each row
                        % corresponding to the assumed node number and
                        % corresponding to the nodes used in the "faces" input
        % Outputs:
                % quality_indexes: A vector of (3M x 1) containing the
                        % angles of each of the sides of every face of a mesh.
                % min_angles: A vector of (M x 1) containing the
                        %  minimum angles of each of every face of a mesh.
        quality_indexes=zeros(size(faces,1)*3,1);
        min_angles=zeros(size(faces,1),1);
        for count_face=1:size(faces,1)
                nodel=faces(count_face,:);
                face_nodes=nodes(nodel,:);
                vec1=face_nodes(2,:)-face_nodes(1,:);
                vec1_norm=vec1/norm(vec1);
                vec2=face_nodes(3,:)-face_nodes(2,:);
                vec2_norm=vec2/norm(vec2);
                vec3=face_nodes(3,:)-face_nodes(1,:);
                vec3_norm=vec3/norm(vec3);
                angle1=rad2deg(acos(dot(vec1_norm,vec3_norm)));
                angle2=rad2deg(acos(dot(vec2_norm,-vec1_norm)));
                angle3=180-angle1-angle2;
                quality_indexes(count_face*3-2)=angle1;
                quality_indexes(count_face*3-1)=angle2;
                quality_indexes(count_face*3)=angle3;
                min_angles(count_face)=min([angle1,angle2,angle3]);
        end

end