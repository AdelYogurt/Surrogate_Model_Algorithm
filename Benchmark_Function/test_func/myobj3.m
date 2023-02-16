function f = myobj3(x)
%-------Ackley Problem---------%

f = -20*exp(-0.2*sqrt(sum(x.^2)/size(x,2))) - exp(sum(cos(2*pi*x))/size(x,2)) +20 + exp(1);

end
