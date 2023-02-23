clc;
% clear;
close all hidden;

sample_number=20;
variable_number=3;
low_bou=-0*ones(1,variable_number);
up_bou=1*ones(1,variable_number);
% x_exist_list=[0.1,0.5,0.3;0.6,0.3,0.6];
% x_exist_list=[0.1,0.5;0.6,0.3];
x_exist_list=[];

radius=1;

A=[];
B=[];
Aeq=[];
Beq=[];

% cheapcon_function=@(x) sum(x.^2)-radius*radius;
% % normalize cheapcon_function
% % optiaml
% ga_option=optimoptions('ga','FunctionTolerance', 1e-2,...
%     'PopulationSize',max(10,2*variable_number),...
%     'MaxGenerations',50,...
%     'display','none');
% [~,fval_max]=ga...
%     (@(x) -cheapcon_function(x),variable_number,A,B,Aeq,Beq,low_bou,up_bou,[],ga_option);
% if fval_max < 0
%     cheapcon_function=@(x) -cheapcon_function(x)/fval_max;
% end

cheapcon_function=[];

tic
[X,X_supply,distance_min_nomlz]=getLatinHypercube...
    (sample_number,variable_number,x_exist_list,...
    low_bou,up_bou,cheapcon_function);
toc
line(X(:,1),X(:,2),X(:,3),'linestyle','none','Marker','o','color','b');
hold on;
% rectangle('Position',[-radius,-radius,2*radius,2*radius],'Curvature',[1 1])

bou=[low_bou(1:2);up_bou(1:2)];
axis(bou(:));
grid on;
hold off;

tic;
x_lhs=lhsdesign(sample_number,variable_number,"iterations",800);
toc;

line(x_lhs(:,1),x_lhs(:,2),x_lhs(:,3),'linestyle','none','Marker','.','color','r');

view(3);

function [X,X_new,distance_min_nomlz]=getLatinHypercube...
    (sample_number,variable_number,X_exist,...
    low_bou,up_bou,cheapcon_function)
% generate sample sequence latin hypercube
% more uniform point distribution by simulating particle motion
% sample number is total point in area
% default low_bou is 0, up_bou is 1, cheapcon_function is []
% low_bou and up_bou is colume vector
% x in x_exist_list, x_list, supply_x_list is row vector
% x_exist_list should meet bou
%
% Copyright 2022 Adel
%
if nargin < 6
    cheapcon_function=[];
    if nargin < 5
        if nargin < 3
            X_exist=[];
            if nargin < 2
                error('getLatinHypercube: lack variable_number');
            end
        end
        low_bou=zeros(1,variable_number);
        up_bou=ones(1,variable_number);
    end
end
iteration_max=100;

% check x_exist_list if meet boundary
if ~isempty(X_exist)
    index=find(X_exist < low_bou);
    index=[index,find(X_exist > up_bou)];
    if ~isempty(index)
        error('getLatinHypercube: x_exist_list range error');
    end
    if size(X_exist,2) ~= variable_number
        error('getLatinHypercube: x_exist_list variable_number error');
    end
    X_exist_nomlz=(X_exist-low_bou)./(up_bou-low_bou);
else
    X_exist_nomlz=[];
end

% check input
if sample_number < 0
    X=[];
    X_new=[];
    distance_min_nomlz=[];
    return;
end

% check x_new_number
x_new_number=sample_number-size(X_exist,1);
if x_new_number < 0
    X=X_exist;
    X_new=[];
    distance_min_nomlz=getMinDistance(X_exist_nomlz);
    return;
end

low_bou_nomlz=zeros(1,variable_number);
up_bou_nomlz=ones(1,variable_number);

% get initial X_new_nomalize by lhsdesign
X_new_nomlz=rand(x_new_number,variable_number);
distance_min_nomlz=getMinDistance([X_new_nomlz;X_exist_nomlz]);

% x is nomalize, so constraint function should change
if ~isempty(cheapcon_function)
    cheapcon_function=@(x) ...
        max(max(sample_number*cheapcon_function(x.*(up_bou-low_bou)+low_bou)+1,0),[],1);
end

% pic_num=1;

iteration=0;
fval_list=zeros(x_new_number,1);
gradient_list=zeros(x_new_number,variable_number);
while iteration < iteration_max
%     scatter(X_new_nomlz(:,1),X_new_nomlz(:,2));
%     bou=[low_bou_nomlz(1:2);up_bou_nomlz(1:2)];
%     axis(bou(:));
%     grid on;
%     
%     radius=1;
%     hold on;
%     rectangle('Position',[-radius,-radius,2*radius,2*radius],'Curvature',[1 1])
%     hold off;
%     
%     drawnow;
%     F=getframe(gcf);
%     I=frame2im(F);
%     [I,map]=rgb2ind(I,256);
%     if pic_num == 1
%         imwrite(I,map,'show_trajectory_constrain.gif','gif','Loopcount',inf,'DelayTime',0.1);
%     else
%         imwrite(I,map,'show_trajectory_constrain.gif','gif','WriteMode','append','DelayTime',0.1);
%     end
%     pic_num = pic_num + 1;
    
    % change each x place by newton methods
    for x_index=1:x_new_number
%         draw_function=@(x) objectFunctionXPlace...
%             (x,[X_new_nomlz(1:x_index-1,:);X_new_nomlz(x_index+1:end,:);X_exist_nomlz],...
%             sample_number,variable_number,low_bou_nomlz-0.1/variable_number,up_bou_nomlz+0.1/variable_number,cheapcon_function);
%         drawFunction(draw_function,low_bou_nomlz,up_bou_nomlz);
        
        % get gradient
        [fval_list(x_index,1),gradient_list(x_index,:)]=objectFunctionXPlace...
            (X_new_nomlz(x_index,:),[X_new_nomlz(1:x_index-1,:);X_new_nomlz(x_index+1:end,:);X_exist_nomlz],...
            sample_number,variable_number,low_bou_nomlz-0.1/variable_number,up_bou_nomlz+0.1/variable_number,cheapcon_function);
        
%         [gradient,hessian]=differ(draw_function,X_new_nomlz(x_index,:));
    end
    
    % normalize fval
    fval_list=fval_list/max(fval_list);
    for x_index=1:x_new_number
        C=fval_list(x_index,1)*distance_min_nomlz*(1-iteration/iteration_max);
        x=X_new_nomlz(x_index,:)+...
            -gradient_list(x_index,:)/...
            norm(gradient_list(x_index,:))*C;
        x=min(x,up_bou_nomlz);
        x=max(x,low_bou_nomlz);
        X_new_nomlz(x_index,:)=x;
    end
    
    iteration=iteration+1;
end
distance_min_nomlz=getMinDistance([X_new_nomlz;X_exist_nomlz]);
X_new=X_new_nomlz.*(up_bou-low_bou)+low_bou;
X=[X_new;X_exist];

    function [fval,gradient]=objectFunctionXPlace...
            (x,X_surplus,sample_number,variable_number,low_bou,up_bou,cheapcon_function)
        % function describe distance between X and X_supply
        % x is colume vector and X_surplus is matrix which is num-1 x var
        % low_bou_limit__ and up_bou_limit__ is colume vector
        % variable in colume
        %
        a__=10/variable_number;
        a_bou__=30/sample_number;
        
        sign__=((x > X_surplus)-0.5)*2;
        
        xi__=-a__*(x-X_surplus).*sign__;
        sum_xi__=sum(xi__,2);
        psi__=a__*(low_bou-x)*a_bou__;
        zeta__=a__*(x-up_bou)*a_bou__;
        
%         exp_xi__=exp(xi__);
        exp_sum_xi__=exp(sum_xi__);
        exp_psi__=exp(psi__);
        exp_zeta__=exp(zeta__);
        
        % get fval
        fval=sum(exp_sum_xi__,1)+...
            sum(exp_psi__+exp_zeta__,2);
        
        % get gradient
        gradient=sum(-a__*sign__.*exp_sum_xi__,1)+...
            -a__*exp_psi__*a_bou__+...
            a__*exp_zeta__*a_bou__;
        
        if ~isempty(cheapcon_function)
            fval_con=cheapcon_function(x);
            fval=fval+fval_con;
            [gradient_con]=differ...
                (cheapcon_function,x,fval_con,variable_number);
            gradient=gradient+gradient_con;
        end
        
        function [gradient]=differ(differ_function,x,fval,variable_number,step)
            % differ function to get gradient and hessian
            % gradient is rank vector
            %
            if nargin < 5
                step=1e-6;
            end
            fval__=zeros(variable_number,2); % backward is 1, forward is 2
            gradient=zeros(1,variable_number);
            
            % fval and gradient
            for variable_index__=1:variable_number
                x_forward__=x;
                x_backward__=x;
                x_backward__(variable_index__)=x_backward__(variable_index__)-step;
                fval__(variable_index__,1)=differ_function(x_backward__);
                
                x_forward__(variable_index__)=x_forward__(variable_index__)+step;
                fval__(variable_index__,2)=differ_function(x_forward__);
                
                gradient(variable_index__)=...
                    (fval__(variable_index__,2)-fval__(variable_index__,1))/2/step;
            end
        end
    end
    function distance_min__=getMinDistance(x_list__)
        % get distance min from x_list
        %
        
        % sort x_supply_list_initial to decrese distance calculate times
        x_list__=sortrows(x_list__,1);
        sample_number__=size(x_list__,1);
        variable_number__=size(x_list__,2);
        distance_min__=variable_number__;
        for x_index__=1:sample_number__
            x_curr__=x_list__(x_index__,:);
            x_next_index__=x_index__ + 1;
            % first dimension only search in min_distance
            search_range__=variable_number__;
            while x_next_index__ <= sample_number__ &&...
                    (x_list__(x_next_index__,1)-x_list__(x_index__,1))^2 ...
                    < search_range__
                x_next__=x_list__(x_next_index__,:);
                distance_temp__=sum((x_next__-x_curr__).^2);
                if distance_temp__ < distance_min__
                    distance_min__ = distance_temp__;
                end
                if distance_temp__ < search_range__
                    search_range__ = distance_temp__;
                end
                x_next_index__=x_next_index__+1;
            end
        end
        distance_min__=sqrt(distance_min__);
    end
end
