function [con,coneq]=functionG06NonlconLow(x)
[con,coneq]=functionG06Nonlcon(x);
con=0.9*con-0.05;
coneq=0.9*coneq-0.05;
end