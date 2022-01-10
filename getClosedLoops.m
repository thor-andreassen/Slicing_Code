function [loop_segments,closed_loops,num_sets]=getClosedLoops(plane_segment_ids)
%% function to calculate the set of closed loops based on an array of edges
        % Created by Thor Andreassen
        % 1/10/22
        %
        % This function takes in a matrix containing all of the edges at a
        % current "slice" as node number connections for the ends of each
        % line segment, and determines the set of closed loops for each set
        % of segments.
        %
        %
        % Inputs:
                % plane_segment_ids: A matrix of (n x 2) integers,
                        % containing the node numbers at the end of each line
                        % segment. With each row corresponding to another line
                        % segment, and each column corresponding to one of the
                        % node numbers on the ends of the line segment
        %
        %
        % Outputs:
                % loop_segments: A cell array containing each of the loops
                        % as separe array matrices in each cell.
                % closed_loops: An array corresponding to the number of
                        % loops in the loop_segments output, where a 1 is present
                        % if the loop is closed, and a 0 if it is not
                        % closed.
                % num_sets: an integer containing the number of closed loop
                        % sets.
        closed_loops=[];       
        plane_segment_ids2=plane_segment_ids;
        loop_segments={};
        current_loop_segments={};
        counter=1;
        max_count=75;
        loop_count=1;
        while ~isempty(plane_segment_ids2) && counter<=max_count
                current_loop_segments={};
                current_node1=plane_segment_ids2(1,1);
                current_node2=plane_segment_ids2(1,2);
                start_node=current_node1;
                [temp_row,~]=find(plane_segment_ids2==current_node2);
                row_id=1;
                current_loop_segments=[current_loop_segments,plane_segment_ids2(row_id,:)];
                for counti=1:length(temp_row)
                        if ~((plane_segment_ids2(temp_row(counti),1)==current_node1 && ...
                                        plane_segment_ids2(temp_row(counti),2)==current_node2) || ...
                                        (plane_segment_ids2(temp_row(counti),1)==current_node2 && ...
                                        plane_segment_ids2(temp_row(counti),2)==current_node1))
                                row_id=counti;
                        end
                end
                row_id=temp_row(row_id);
                current_node1=current_node2;
                if plane_segment_ids2(row_id,1)==current_node1
                        next_node2=plane_segment_ids2(row_id,2);
                else
                        next_node2=plane_segment_ids2(row_id,1);
                end
                current_loop_segments=[current_loop_segments,plane_segment_ids2(row_id,:)];
                current_node2=next_node2;
                
                inner_loop_counter=1;
                while (current_node2~=start_node &&...
                                ~isempty(find(plane_segment_ids2==current_node2, 1))) &&...
                                inner_loop_counter<size(plane_segment_ids2,1)
                        [temp_row,~]=find(plane_segment_ids2==current_node2);
                        row_id=1;
                        for counti=1:length(temp_row)
                                if ~((plane_segment_ids2(temp_row(counti),1)==current_node1 && ...
                                                plane_segment_ids2(temp_row(counti),2)==current_node2) || ...
                                                (plane_segment_ids2(temp_row(counti),1)==current_node2 && ...
                                                plane_segment_ids2(temp_row(counti),2)==current_node1))
                                        row_id=counti;
                                end
                        end

                        row_id=temp_row(row_id);
                        current_node1=current_node2;
                        if plane_segment_ids2(row_id,1)==current_node1
                                next_node2=plane_segment_ids2(row_id,2);
                        else
                                next_node2=plane_segment_ids2(row_id,1);
                        end
                        current_loop_segments=[current_loop_segments,plane_segment_ids2(row_id,:)];
                        current_node2=next_node2;
                        inner_loop_counter=inner_loop_counter+1;
                end
                if (current_node2==start_node)
                        closed_loops=[closed_loops,1];
                        temp_mat=unique(cell2mat([current_loop_segments]),'stable');
                        temp_mat=[temp_mat,temp_mat(1)];
                else
                        closed_loops=[closed_loops,0];
                        temp_mat=unique(cell2mat([current_loop_segments]),'stable');
                end
                
                loop_segments{loop_count}=temp_mat;
                loop_count=loop_count+1;
                for count_segments=1:length(current_loop_segments)
                        [common_id]=find(plane_segment_ids2==current_loop_segments{count_segments});
                        if ~isempty(common_id)
                                if size(plane_segment_ids2,1)>=common_id(1)
                                        plane_segment_ids2(common_id(1),:)=[];
                                end
                        end
                end
                counter=counter+1;
        end
        num_sets=counter-1;
end