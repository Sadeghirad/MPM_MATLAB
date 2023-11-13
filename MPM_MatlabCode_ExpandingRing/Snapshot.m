
hold on
grid off

% draw grid cells
for i_g1 = 1:numelem

    patch(x(1,conn(:,i_g1)),x(2,conn(:,i_g1)),'-','LineWidth',0.2, ...
        'FaceColor','none','EdgeColor','k')

end

% draw particle domains
for i_g1 = 1:numpar

    if strcmp(interpolator,'BSMPM')

        ParVec1 = InitParVec1(:,i_g1);
        ParVec2 = InitParVec2(:,i_g1);

    elseif strcmp(interpolator,'CPDI')  || strcmp(interpolator,'BSCPDI')

        ParVec1 = F(:,:,i_g1)*InitParVec1(:,i_g1);
        ParVec2 = F(:,:,i_g1)*InitParVec2(:,i_g1);

    end

    xEx = zeros(2,4);
    xEx(:,1) = [xp(1,i_g1)-ParVec1(1)-ParVec2(1);xp(2,i_g1)-ParVec1(2)-ParVec2(2)];
    xEx(:,2) = [xp(1,i_g1)+ParVec1(1)-ParVec2(1);xp(2,i_g1)+ParVec1(2)-ParVec2(2)];
    xEx(:,3) = [xp(1,i_g1)+ParVec1(1)+ParVec2(1);xp(2,i_g1)+ParVec1(2)+ParVec2(2)];
    xEx(:,4) = [xp(1,i_g1)-ParVec1(1)+ParVec2(1);xp(2,i_g1)-ParVec1(2)+ParVec2(2)];

    patch(xEx(1,:),xEx(2,:),'-','LineWidth',0.5,'FaceColor','c','EdgeColor','k')

end

% draw particles
plot(xp(1,:),xp(2,:),'bo','LineWidth',0.2,'MarkerEdgeColor','b', ...
    'MarkerFaceColor','b','MarkerSize',2)
axis equal
axis off