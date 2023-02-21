classdef BenchmarkFunction < handle
    methods % unconstraint single objective function
        function fval=single2DObject(self,x)
            % Rosenbrock problem
            %
            % variable_number=2;
            % object_function=@(x) benchmark.single2DObject(x);
            % object_function_LF=@(x) benchmark.single2DObjectLow(x);
            % A=[];
            % B=[];
            % Aeq=[];
            % Beq=[];
            % low_bou=[-5,-5];
            % up_bou=[5,5];
            % nonlcon_function=[];
            % nonlcon_function_LF=[];
            % cheapcon_function=[];
            %
            % x_best=[-3.6483,-0.0685] fval_best=272.5563
            %
            x1=x(1);
            x2=x(2);
            fval=(x1^2+x2-11)^2+(x2^2+x1+20)^2;
        end
        function fval=single2DObjectLow(self,x)
            x1=x(1);
            x2=x(2);
            fval=self.single2DObject([0.9*x1;0.8*x2])-(x1+1)^2;
        end
        function fval=singleA10Object(self,x)
            x=x(:);
            s=[1.3;0.1;1.4;0.8;1.7;1;1.5;0.6;2;0.4];
            fval=-20*exp(-0.2*sqrt(sum((x-s).^2)/10))-exp(sum(cos(2*1.3*pi*(x-s))/10))+20+exp(1);
        end
        function fval=singleA10ObjectLow(self,x)
            fval=-20*exp(-0.2*sqrt(sum(x.^2)/10))-exp(sum(cos(2*pi*x)/10))+20+exp(1);
        end
        function fval=singleA20Object(self,x)
            s=[1.3;0.1;1.4;0.8;1.7;1;1.5;0.6;2;0.4;1.3;0.3;1.5;0.9;1.9;1.1;1.7;0.7;2.1;0.5];
            fval=-20*exp(-0.2*sqrt(sum((x-s).^2)/10))-exp(sum(cos(2*1.3*pi*(x-s))/10))+20+exp(1);
        end
        function fval=singleA20ObjectLow(self,x)
            fval=-20*exp(-0.2*sqrt(sum(x.^2)/10))-exp(sum(cos(2*pi*x)/10))+20+exp(1);
        end
        function fval=singleBRObject(self,x)
            % Branin function
            % variable_number=2;
            % object_function=@(x) functionBRObject(x);
            % A=[];
            % B=[];
            % Aeq=[];
            % Beq=[];
            % low_bou=[-5;10];
            % up_bou=[0;15];
            % nonlcon_function=[];
            % cheapcon_function=[];
            % model_function=[];
            %
            % x_best=[-3.1416;12.2750] fval_best=0.3979
            %
            x1=x(1);
            x2=x(2);
            fval=(x2-5.1/4/pi/pi*x1^2+5/pi*x1-6)^2+10*(1-1/8/pi)*cos(x1)+10;
        end
        function fval=singleCOLObject(self,x)

            x1=x(1);
            x2=x(2);
            x3=x(3);
            x4=x(4);
            fval=100*(x1^2-x2)^2+(x1-1)^2+(x3-1)^2+90*(x3^2-x4)^2+...
                10.1*((x2-1)^2-(x4-1)^2)+19.8*(x2-1)*(x4-1);
        end
        function fval=singleCOLObjectLow(self,x)
            fval=func_H2_COL_H([0.8;0.8;0.5;0.5].*x);
        end
        function fval=singleCOSObject(self,x)
            x=x(:);
            par=[2,-2];
            fval=sum(sum((x-par').^2-10*cos(2*pi*x)+10));
        end
        function fval=singleDP20Object(self,x)
            s=[1.8;0.5;2;1.2;0.4;0.2;1.4;0.3;1.6;0.6;0.8;1;1.3;1.9;0.7;1.6;0.3;1.1;2;1.4];
            fval=func_G2_DP20_L(x-s);
        end
        function fval=singleDP20ObjectLow(self,x)
            fval=(x(1)-1)^2+sum(linspace(2,20,19)'.*(2*x(2:end).^2-x(1:end-1)).^2);
        end
        function fval=singleE20Object(self,x)
            ss=[1.8;0.4;2;1.2;1.4;0.6;1.6;0.2;0.8;1;1.3;1.1;2;1.4;0.5;0.3;1.6;0.7;0.3;1.9];
            sh=[0.3;0.4;0.2;0.6;1;0.9;0.2;0.8;0.5;0.7;0.4;0.3;0.7;1;0.9;0.6;0.2;0.8;0.2;0.5];
            fval=sum(linspace(1,20,20)'.*sh.*(x-ss).^2);
        end
        function fval=singleE20ObjectLow(self,x)
            fval=sum(linspace(1,20,20)'.*x.^2);
        end
        function fval=singleG01Object(self,x)
            % object_function=@benchmark.singleG01Object;
            % object_function_low=@benchmark.singleG01ObjectLow;
            % A=[ 2   2   0   0   0   0   0   0   0   1   1   0   0;
            %     2   0   2   0   0   0   0   0   0   1   0   1  0;
            %     0   2   2   0   0   0   0   0   0   0   1   1  0;
            %     -8  0   0   0   0   0   0   0   0   1   0   0   0;
            %     0   -8  0   0   0   0   0   0   0   0   1   0   0;
            %     0   0   0   -2  -1  0   0   0   0   1   0   0   0;
            %     0   0   0   0   0   -2  -1  0   0   0   1   0   0;
            %     0   0   0   0   0   0   0   -2  -1  0   0   1   0;
            %     ];
            % B=[10;10;10;0;0;0;0;0];
            % Aeq=[];
            % Beq=[];
            % low_bou=zeros(13,1);
            % up_bou=ones(13,1);
            % up_bou(10:12)=100;
            % nonlcon_function=[];
            %
            sigma1=0;
            for i=1:4
                sigma1=sigma1+x(i);
            end
            sigma2=0;
            for i=1:4
                sigma2=sigma2+x(i)^2;
            end
            sigma3=0;
            for i=5:13
                sigma3=x(i)+sigma3;
            end
            fval=5*sigma1-5*sigma2-sigma3;
        end
        function fval=singleG01ObjectLow(self,x)
            fval=self.singleG01Object(x)*0.9+0.5;
        end
        function fval=singleForresterObject(self,x)
            fval=((x.*6-2).^2).*sin(x.*12-4);
        end
        function fval=singleForresterObjectLow(self,x)
            A=0.5;B=10;C=-5;
            fval=A*self.singleForresterObject(x)+B*(x-0.5)+C;
        end
        function fval=singleGFObject(self,x)
            x=x(:);
            c1=1.5;
            c2=2.25;
            c3=2.625;
            x1=x(1);
            x2=x(2);
            u1=c1-x1*(1-x2^1);
            u2=c2-x1*(1-x2^2);
            u3=c3-x1*(1-x2^3);
            fval=u1^2+u2^2+u3^2;
        end
        function fval=singleGPObject(self,x)
            % variable_number=2;
            % object_function=@(x) benchmark.singleGPObject(x);
            % A=[];
            % B=[];
            % Aeq=[];
            % Beq=[];
            % low_bou=[-2;-2];
            % up_bou=[2;2];
            % nonlcon_function=[];
            % cheapcon_function=[];
            %
            % x_best=[0;-1] fval_best=3
            %
            if (size(x,2) ~= 2)
                error('singleGPObject: variable must be 2');
            end
            x1=x(:,1);
            x2=x(:,2);
            fval=(1+(x1+x2+1).^2.*...
                (19-14*x1+3*x1.^2-14*x2+6*x1.*x2+3*x2.^2)).*...
                (30+(2*x1-3*x2).^2.*(18-32*x1+12*x1.^2+48*x2-36*x1.*x2+27*x2.^2));
        end
        function fval=singleHNObject(self,x)
            % Hartman function
            %
            % object_function=@(x) benchmark.singleHNObject(x);
            % variable_number=6;
            % A=[];
            % B=[];
            % Aeq=[];
            % Beq=[];
            % low_bou=zeros(1,variable_number);
            % up_bou=ones(1,variable_number);
            % nonlcon_function=[];
            % cheapcon_function=[];
            % model_function=[];
            %
            % x_best=[0.2017;0.1500;0.4769;0.2753;0.3117;0.6573] fval_best=-3.3224
            %
            x=x(:)';
            coe=[1	10	3	17	3.5	1.7	8	1	0.1312	0.1696	0.5569	0.0124	0.8283	0.5886;
                2	0.05	10	17	0.1	8	14	1.2	0.2329	0.4135	0.8307	0.3736	0.1004	0.9991;
                3	3	3.5	1.7	10	17	8	3	0.2348	0.1451	0.3522	0.2883	0.3047	0.6650;
                4	17	8	0.05	10	0.1	14	3.2	0.4047	0.8828	0.8732	0.5743	0.1091	0.0381;];

            alpha=coe(:,2:7);
            c=coe(:,8);
            p=coe(:,9:14);

            fval=0;
            for i=1:4
                hari=0;
                for j=1:6
                    hari=alpha(i,j).*(x(:,j)-p(i,j)).^2+hari;
                end
                fval=c(i)*exp(-hari)+fval;
            end
            fval=-fval;
        end
        function fval=singlePKObject(self,x)
            % variable_number=2;
            % object_function=@(x) benchmark.singlePKObject(x);
            % A=[];
            % B=[];
            % Aeq=[];
            % Beq=[];
            % low_bou=[-3,-3];
            % up_bou=[3,3];
            % nonlcon_function=[];
            % cheapcon_function=[];
            % model_function=[];
            %
            % x_best=[0.2283,-1.6255], fval_best=-6.5511
            %
            x1=x(1);
            x2=x(2);
            fval=3*(1-x1).^2.*exp(-(x1.^2) - (x2+1).^2) ...
                - 10*(x1/5 - x1.^3 - x2.^5).*exp(-x1.^2-x2.^2) ...
                - 1/3*exp(-(x1+1).^2 - x2.^2) ;
        end
        function fval=singlePloyObject(self,x)
            fval=520*x.^5-1000*x.^4+585*x.^3-79*x.^2-20*x+4;
        end
        function fval=singleR10Object(self,x)
            fval=sum(100*(x(2:10)-x(1:9).^2).^2+(x(1:9)-1).^2);
        end
        function fval=singleROSObject(self,x)
            fval=sum((100*(x(2:4)-x(1:3).^2).^2-(x(1:3)-1).^2).^2);
        end
        function fval=singleROSObjectLow(self,x)
            fval=sum((100*(0.5*x(2:4)-0.6*x(1:3).^2).^2-(0.5*x(1:3)-0.5).^2).^2);
        end
        function fval=singleRSObject(self,x)
            x1=x(1);
            x2=x(2);
            fval=x1^2+x2^2-cos(18*x1)-cos(18*x2);
        end
        function fval=singleSCObject(self,x)
            x1=x(1);
            x2=x(2);
            fval=4*x1^2-2.1*x1^4+x1^6/3+x1*x2-4*x2^2+4*x2^4;
        end
        function fval=singleST5Object(self,x)
            s=[0.28;0.59;0.47;0.16;0.32];
            fval=func_H2_ST5_L(x-s);
        end
        function fval=singleST5ObjectLow(self,x)
            fval=0.5*sum(x.^4-16*x.^2+5*x);
        end

    end
    methods % constraint single objective function
        function [con,coneq]=singleICENonlcon(self,x)
            b=x(1);
            c_r=x(2);
            d_E=x(3);
            d_I=x(4);
            omega=x(5)*1e-3;

            K_1=1.2;
            K_2=2;
            K_3=0.82;
            K_4=(-1e-12+30.99)/37.34; K_5=0.89;
            K_6=0.6;
            K_7=6.5;
            K_8=230.5;
            L_1=400;
            L_2=200;
            rou=1.225;
            gama=1.33;
            V=1.859*1e6;
            Q=43958;
            N_c=4;
            C_s=0.44;
            A_f=14.6;

            S_V=0.83*((8+4*c_r)+1.5*(c_r-1)*(pi*N_c/V)*b^3)/((2+c_r)*b);
            eta_tw=0.8595*(1-c_r^(-0.33))-S_V;

            g1=K_1*N_c*b-L_1;
            g2=(4*K_2*V/pi/N_c/L_2)^0.5-b;
            g3=d_I+d_E-K_3*b;
            g4=K_4*d_I-d_E;
            g5=d_E-K_5*d_I;
            g6=9.428*1e-5*(4*V/pi/N_c)*(omega/d_I^2)-K_6*C_s;
            g7=c_r-13.2+0.045*b;
            g8=omega-K_7;
            g9=3.6*1e6-K_8*Q*eta_tw;

            con=[g1;g2;g3;g4;g5;g6;g7;g8;g9];
            coneq=[];
        end
        function fval=singleICEObject(self,x)

            b=x(1);
            c_r=x(2);
            d_E=x(3);
            d_I=x(4);
            omega=x(5)*1e-3;

            K_0=1/120;
            rou=1.225;
            gama=1.33;
            V=1.859*1e6;
            Q=43958;
            N_c=4;
            C_s=0.44;
            A_f=14.6;

            if omega > 5.25
                eta_vb=1.067-0.038*exp(omega-5.25);
            else
                eta_vb=0.637+0.13*omega-0.014*omega^2+0.00066*omega^3;
            end
            eta_tad=0.8595*(1-c_r^(-0.33));
            eta_V=eta_vb*(1+5.96*1e-3*omega^2)/(1+((9.428*1e-5)*(4*V/pi/N_c/C_s)*(omega/d_I^2))^2);
            S_V=0.83*((8+4*c_r)+1.5*(c_r-1)*(pi*N_c/V)*b^3)/((2+c_r)*b);
            eta_t=eta_tad-S_V*(1.5/omega)^0.5;
            V_P=(8*V/pi/N_c)*omega*b^(-2);
            FMEP=4.826*(c_r-9.2)+(7.97+0.253*V_P+9.7*(1e-6)*V_P^2);

            fval=K_0*(FMEP-(rou*Q/A_f)*eta_t*eta_V)*omega;
        end

        function [con,coneq]=singlePVDNonlcon(self,x)
            x1=x(1);x2=x(2);x3=x(3);x4=x(4);
            g3=-pi*x3^2-4/3*pi*x3^3+1296000;
            g4=x4-240;
            g5=1.1-x1;
            g6=0.6-x2;
            con=[g3;g4;g5;g6];
            coneq=[];
        end
        function fval=singlePVDObject(self,x)
            x1=x(1);x2=x(2);x3=x(3);x4=x(4);
            fval=0.6224*x1*x3*x4+1.7781*x2*x3^2+3.1661*x1^2*x4+19.84*x1^2*x3;
        end

        function [con,coneq]=singlePVD4Nonlcon(self,x)
            x1=x(1);x2=x(2);x3=x(3);x4=x(4);
            g3=-pi*x3^2*x4-4/3*pi*x3^3+1296000;
            if g3 >=0
                g3=log(1+g3);
            else
                g3=-log(1-g3);
            end
            con=g3;
            coneq=[];
        end
        function [con,coneq]=singlePVD4NonlconLow(self,x)
            [con,coneq]=singlePVDNonlcon(self,x);
            con=0.9*con-0.05;
            coneq=0.9*coneq-0.05;
        end
        function fval=singlePVD4Object(self,x)
            x1=x(1);x2=x(2);x3=x(3);x4=x(4);
            fval=0.6224*x1*x3*x4+1.7781*x2*x3^2+3.1661*x1^2*x4+19.84*x1^2*x3;
        end
        function fval=singlePVD4ObjectLow(self,x)
            fval=singlePVDObject(self,x)*0.9+0.5;
        end

        function [con,coneq]=singleG06Nonlcon(self,x)
            g1=-(x(1)-5)^2-(x(2)-5)^2+100;
            g2=(x(1)-6)^2+(x(2)-5)^2-82.81;
            con=[g1;g2];
            coneq=[];
        end
        function [con,coneq]=singleG06NonlconLow(self,x)
            [con,coneq]=singleG06Nonlcon(self,x);
            con=0.9*con-0.05;
            coneq=0.9*coneq-0.05;
        end
        function fval=singleG06Object(self,x)
            fval=(x(1)-10)^3+(x(2)-20)^3;
        end
        function fval=singleG06ObjectLow(self,x)
            fval=0.9*singleG06Object(self,x)+0.5;
        end
    end
    methods % unconstraint mulit objective function
        function fval=multiZDT1Object(self,x)
            variable_number=10;
            fval=zeros(2,1);
            fval(1) = x(1);
            g=1+9*(sum(x(2:variable_number))/(variable_number-1));
            fval(2) = g * (1 - sqrt(x(1) / g));
        end
        function fval=multiZDT2Object(self,x)
            variable_number=10;
            fval=zeros(2,1);
            fval(1)=x(1);
            g=1+9*(sum(x(2:variable_number))/(variable_number-1));
            fval(2)=g*(1-(x(1)/g)^2);
        end
        function fval=multiZDT3Object(self,x)
            variable_number=10;
            fval=zeros(2,1);
            fval(1)=x(1);
            g=1+9*(sum(x(2:variable_number))/(variable_number-1));
            fval(2)=g*(1-sqrt(x(1)/g)-(x(1)/g)*sin(10*pi*x(1)));
        end
    end
    methods
        function [con,coneq]=multiTNKNonlcon(self,x)
            con=zeros(2,1);
            x1=x(1);
            x2=x(2);
            if x1 == 0 && x2 == 0
                con(1)=-(x1^2+x2^2-1-0.1);
            else
                con(1)=-(x1^2+x2^2-1-0.1*cos(16*atan(x1/x2)));
            end
            con(2)=(x1-0.5)^2+(x2-0.5)^2-0.5;
            coneq=[];
        end
        function fval=multiTNKObject(self,x)
            % TNK problem
            % variable_number is 2
            %
            % object_function=@(x) benchmark.multiTNKObject(self,x);
            % variable_number=2;
            % low_bou=zeros(1,2);
            % up_bou=ones(1,2)*pi;
            % nonlcon_function=@(x) multiTNKNonlcon(self,x);
            % cheapcon_function=[];
            %
            fval=zeros(2,1);
            fval(1)=x(1);
            fval(2)=x(2);
        end
    end
end