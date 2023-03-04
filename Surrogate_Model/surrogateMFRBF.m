clc;
clear;
close all hidden;

% load('Forrester.mat')
% [predict_function_HK, HK_model] = interpHieraKrigingPreModel...
%     (XHF, YHF, XLF, YLF);
% [predict_function, K_model] = interpKrigingPreModel(XHF, YHF);
% [y_HK] = predict_function_HK(x);
% [y] = predict_function(x);
% line(x, y_real, 'Color', 'k', 'LineStyle', '-')
% line(x, y_real_low, 'Color', 'k', 'LineStyle', '--');
% line(x, y_HK, 'Color', 'b', 'LineStyle', '-')
% line(x, y, 'Color', 'b', 'LineStyle', '--');
% legend({'singleForresterObject', 'singleForresterObjectLow', 'HK', 'K'})

load('2DMF.mat');
[predict_function_HK, MFRBF_model] = interpMultiRadialBasisPreModel(XHF, YHF, [], XLF, YLF, []);
interpVisualize(MFRBF_model, low_bou, up_bou);

function [predict_function_MFRBF, MFRBF_model] = interpMultiRadialBasisPreModel...
    (XHF, YHF, varargin)
% multi fildelity radial basis function interp pre model function
% XHF, YHF are x_HF_number x variable_number matrix
% XLF, YLF are x_LF_number x variable_number matrix
% aver_X, stdD_X is 1 x x_HF_number matrix
%
% input:
% XHF, YHF, basis_func_HF(can be []), XLF, YLF, basis_func_LF(can be [])
% XHF, YHF, basis_func_HF(can be []), LF_model
%
% output:
% predict_function, HK_model
%
% reference: [1] LIU Y, WANG S, ZHOU Q, et al. Modified Multifidelity
% Surrogate Model Based on Radial Basis Function with Adaptive Scale Factor
% [J]. Chinese Journal of Mechanical Engineering, 2022, 35(1): 77.
%
% Copyright 2023 Adel
%
[x_HF_number, variable_number] = size(XHF);
switch nargin
    case 4
        basis_func_HF = varargin{1};
        LF_model = varargin{2};

        % check whether LF model exist predict_function
        if ~isfield(LF_model, 'predict_function')
            error('interpHieraKrigingPreModel: low fidelity lack predict function');
        end
    case 6
        basis_func_HF = varargin{1};
        XLF = varargin{2};
        YLF = varargin{3};
        basis_func_LF = varargin{4};

        [x_LF_number, variable_number] = size(XLF);

        % first step
        % construct low fidelity model

        % normalize data
        aver_X = mean(XLF);
        stdD_X = std(XLF);
        aver_Y = mean(YLF);
        stdD_Y = std(YLF);
        index__ = find(stdD_X == 0);
        if  ~isempty(index__), stdD_X(index__) = 1; end
        index__ = find(stdD_Y == 0);
        if  ~isempty(index__), stdD_Y(index__) = 1; end

        XLF_nomlz = (XLF-aver_X)./stdD_X;
        YLF_nomlz = (YLF-aver_Y)./stdD_Y;

        if isempty(basis_func_LF)
            c = (prod(max(XLF_nomlz) - min(XLF_nomlz))/x_LF_number)^(1/variable_number);
            basis_func_LF = @(r) exp(-(r.^2)/c);
        end

        % initialization distance of XLF_nomlz
        XLF_dis = zeros(x_LF_number, x_LF_number);
        for variable_index = 1:variable_number
            XLF_dis = XLF_dis + ...
                (XLF_nomlz(:, variable_index) - XLF_nomlz(:, variable_index)').^2;
        end
        XLF_dis = sqrt(XLF_dis);

        [beta_LF, rdibas_matrix_LF] = interpRadialBasis...
            (XLF_dis, YLF_nomlz, basis_func_LF, x_LF_number);

        % initialization predict function
        predict_function_LF = @(X_predict) interpRadialBasisPredictor...
            (X_predict, XLF_nomlz, aver_X, stdD_X, aver_Y, stdD_Y, ...
            x_LF_number, variable_number, beta_LF, basis_func_LF);

        LF_model.X = XLF;
        LF_model.Y = YLF;
        LF_model.radialbasis_matrix = rdibas_matrix_LF;
        LF_model.beta = beta_LF;

        LF_model.aver_X = aver_X;
        LF_model.stdD_X = stdD_X;
        LF_model.aver_Y = aver_Y;
        LF_model.stdD_Y = stdD_Y;
        LF_model.basis_function = basis_func_LF;

        LF_model.predict_function = predict_function_LF;
    otherwise
        error('interpHieraKrigingPreModel: error input');
end
MFRBF_model.LF_model = LF_model;
predict_function_LF = LF_model.predict_function;

% predict LF value at XHF point
YHF_pred = predict_function_LF(XHF);

% nomalizae
YHF_pred_nomlz = (YHF_pred - aver_Y)./stdD_Y;

% second step
% construct MFRBF model

% normalize data
aver_X = mean(XHF);
stdD_X = std(XHF);
aver_Y = mean(YHF);
stdD_Y = std(YHF);
index__ = find(stdD_X == 0);
if ~isempty(index__), stdD_X(index__) = 1;end
index__ = find(stdD_Y == 0);
if ~isempty(index__), stdD_Y(index__) = 1;end
XHF_nomlz = (XHF - aver_X)./stdD_X;
YHF_nomlz = (YHF - aver_Y)./stdD_Y;

if isempty(basis_func_HF)
    c = (prod(max(XHF_nomlz) - min(XHF_nomlz))/x_HF_number)^(1/variable_number);
    basis_func_HF = @(r) exp(-(r.^2)/c);
end

% initialization distance of XHF_nomlz
XHF_dis = zeros(x_HF_number, x_HF_number);
for variable_index = 1:variable_number
    XHF_dis = XHF_dis + ...
        (XHF_nomlz(:, variable_index) - XHF_nomlz(:, variable_index)').^2;
end
XHF_dis = sqrt(XHF_dis);

[beta_HF, rdibas_matrix_HF] = interpMultiRadialBasis...
    (XHF_dis, YHF_nomlz, basis_func_HF, x_HF_number, YHF_pred_nomlz);

% initialization predict function
predict_function_MFRBF = @(X_predict) interpMultiRadialBasisPredictor...
    (X_predict, XHF_nomlz, aver_X, stdD_X, aver_Y, stdD_Y, ...
    x_HF_number, variable_number, beta_HF, basis_func_HF, predict_function_LF);

MFRBF_model.X = XHF;
MFRBF_model.Y = YHF;
MFRBF_model.radialbasis_matrix = rdibas_matrix_HF;
MFRBF_model.beta = beta_HF;

MFRBF_model.aver_X = aver_X;
MFRBF_model.stdD_X = stdD_X;
MFRBF_model.aver_Y = aver_Y;
MFRBF_model.stdD_Y = stdD_Y;
MFRBF_model.basis_function = basis_func_HF;

MFRBF_model.predict_function = predict_function_MFRBF;

% abbreviation:
% num: number, pred: predict, vari: variable
    function [beta, rdibas_matrix] = interpRadialBasis...
            (X_dis, Y, basis_function, x_number)
        % interp polynomial responed surface core function
        % calculation beta
        %
        % Copyright 2022 Adel
        %
        rdibas_matrix = basis_function(X_dis);

        % stabilize matrix
        rdibas_matrix = rdibas_matrix + eye(x_number)*1e-6;

        % solve beta
        beta = rdibas_matrix\Y;
    end

    function [Y_pred] = interpRadialBasisPredictor...
            (X_pred, X_nomlz, aver_X, stdD_X, aver_Y, stdD_Y, ...
            x_num, vari_num, beta, basis_function)
        % radial basis function interpolation predict function
        %
        [x_pred_num, ~] = size(X_pred);

        % normalize data
        X_pred_nomlz = (X_pred - aver_X)./stdD_X;

        % calculate distance
        X_dis_pred = zeros(x_pred_num, x_num);
        for vari_index = 1:vari_num
            X_dis_pred = X_dis_pred + ...
                (X_pred_nomlz(:, vari_index) - X_nomlz(:, vari_index)').^2;
        end
        X_dis_pred = sqrt(X_dis_pred);

        % predict variance
        Y_pred = basis_function(X_dis_pred)*beta;

        % normalize data
        Y_pred = Y_pred*stdD_Y + aver_Y;
    end

    function [beta, rdibas_matrix] = interpMultiRadialBasis...
            (X_dis, Y, basis_function, x_number, YHF_pred)
        % interp polynomial responed surface core function
        % calculation beta
        %
        % Copyright 2022 Adel
        %
        rdibas_matrix = basis_function(X_dis);

        % stabilize matrix
        rdibas_matrix = rdibas_matrix + eye(x_number)*1e-6;

        % add low fildelity value
        H = [rdibas_matrix.*YHF_pred, rdibas_matrix];

        % solve beta
        beta = H'*((H*H')\Y);
    end

    function [Y_pred] = interpMultiRadialBasisPredictor...
            (X_pred, X_nomlz, aver_X, stdD_X, aver_Y, stdD_Y, ...
            x_num, vari_num, beta, basis_function, predict_function_LF)
        % radial basis function interpolation predict function
        %
        [x_pred_num, ~] = size(X_pred);

        % normalize data
        X_pred_nomlz = (X_pred - aver_X)./stdD_X;

        % calculate distance
        X_dis_pred = zeros(x_pred_num, x_num);
        for vari_index = 1:vari_num
            X_dis_pred = X_dis_pred + ...
                (X_pred_nomlz(:, vari_index) - X_nomlz(:, vari_index)').^2;
        end
        X_dis_pred = sqrt(X_dis_pred);

        % predict low fildelity value
        Y_pred_LF = predict_function_LF(X_pred);

        % nomalizae
        Y_pred_LF_nomlz = (Y_pred_LF - aver_Y)./stdD_Y;

        % combine two matrix
        rdibas_matrix_pred = basis_function(X_dis_pred);
        H = [rdibas_matrix_pred.*Y_pred_LF_nomlz, rdibas_matrix_pred];

        % predict variance
        Y_pred = H*beta;

        % normalize data
        Y_pred = Y_pred*stdD_Y + aver_Y;
    end

end
