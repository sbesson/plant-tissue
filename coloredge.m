% Attribute a color to a number ()

function color = coloredge(i)

map = repmat([0 0 0],32,1);
map(1,:) = [1 0 0];
map(2,:) = [0 1 0];
map(3,:) = [0 0 1];
map(4,:) = [1 1 0];
map(5,:) = [0 1 1];
map(6,:) = [1 0 1];
map(7,:) = [1 0.62 .4];
map(8,:) = [.49 1 .83];
map(9,:) = [.5 0 1];

color = map(round(i),:);

% switch round(i)
%     case 0
%         % White
%         color = [0 0 0];
%     case 1
%         %
%         color = [1 0 0];
%     case 2
%         color = [0 0 1];
%     case 3
%         color = [1 1 0];
%     case 4
%         color = [0 1 1];
%     otherwise
%         color = [0 0 0];
% end

end