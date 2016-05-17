BST=10;
Duration=28; % as obtained from NbLambda_SpeedLassoGLM.fig for a single acoustic model
NbLambdas=1:25;% 10 lambdas sounds reasonnable
NbGLMModel=3;
NbWin=length(20:10:150);
NbCells=404;
Time=(BST.*Duration.*NbLambdas.*NbGLMModel + BST.*Duration.*NbGLMModel).*NbWin.*NbCells./(60.*60.*24.*20);
Time=(BST.*Duration.*NbLambdas.*NbGLMModel + BST.*Duration.*NbGLMModel).*NbWin./(60.*60);

% Find how long should be iteration of gradient descent in penalizedwls at
% the maximum
IterationDur = [0.001 0.01 0.1 1 10]; %(in sec)
NbIterationsglmIRLS = 100;
NbLambdas = 10;
NbGLMModel=3;
NbWin=length(20:10:150);
NbCells=404;
Time = (IterationDur.* NbIterationsglmIRLS .* NbLambdas .* NbGLMModel .* BST + IterationDur.* NbIterationsglmIRLS .* NbGLMModel .* BST).*NbWin./(60*60*24);

