function str = format_p(pval)
    if pval < 0.001
        str = '<0.001';
    else
        str = num2str(round(pval, 3));
    end
end