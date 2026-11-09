function value=china_mixed_statistic(T,cortex)
T=prepare_china_table(T);fit=fit_china_mixed(T,cortex);value=fit.beta(:);
end
