function [a,b,c] = numerical_solution(X,Xp,f,rj,beta_post,rout,model)

switch model
    case 'perceptron'
        [a,b] = perceptron_learning(X,Xp,f,rj,beta_post,rout);
        c = 1;
    case 'fmincon'
        [a,b,c] = fmincon_solution(X,Xp,f,rj,beta_post,rout);
    otherwise
        disp('model is not defined')
end



