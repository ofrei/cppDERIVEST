addpath('DERIVESTsuite')

for func_index=1:2
    if func_index==1, f = @exp; r = @(DO)exp(1); x0=1; cfunc='exp'; end
    if func_index==2, f = @(x)(x+power(x,2)+power(x,3)+power(x,4)+power(x,5)); cfunc='poly'; r = @(DO)factorial(DO); x0=0; end;
    StepRatio = 2.0000001;
    Styles = {'Central', 'Forward', 'Backward'};
    for style_index = 1:3
        Style = Styles{style_index};
        for MethodOrder=1:4
            for DerivativeOrder=1:4
                for RombergTerms = 0:3
                    if (lower(Style(1))=='c') && (mod(MethodOrder,2)==1), continue; end  % central methods must be of order 2 or 4
                    if (lower(Style(1))~='c' && MethodOrder==4 && DerivativeOrder == 4 && RombergTerms==3), continue; end  % this case crashes in original DERIVESTsuite implementation
                    [der2, errest2, finaldelta2] = derivest(f, x0, 'DerivativeOrder', DerivativeOrder, 'MethodOrder', MethodOrder, 'RombergTerms', RombergTerms, 'Style', lower(Style));
                    [der,errest,details] = derivest_of(f, x0, 'DerivativeOrder', DerivativeOrder, 'MethodOrder', MethodOrder, 'RombergTerms', RombergTerms, 'Style', lower(Style));
                    assert(der==der2 && errest==errest2 && details.finaldelta==finaldelta2);
                    fprintf('derivest(%s, %i, %i, %i, DerivestStyle_%s, %i, &der, &err, &finaldelta);  ', cfunc, x0, DerivativeOrder, MethodOrder, Style, RombergTerms);
                    fprintf('check(der, err, %.20f, %.20f, %.20f);\n', der, errest, r(DerivativeOrder));
                end
            end
        end
    end
end
%details.der_init(abs(details.der_init)>1e5) = nan;
%min(abs(details.der_init-r))
%min(abs(details.der_romb-r))
%details.errors_romb
%details.finaldelta
%der-1;
%errest;

%[der, errest, details]  = derivest_of(@exp, 1, 'DerivativeOrder', 1, 'MethodOrder', 4, 'RombergTerms', 2, 'Style', 'central')
%derivest_of(@exp, 1, 'DerivativeOrder', 1, 'MethodOrder', 1, 'RombergTerms', 0, 'Style', 'forward')

%    [der,errest,details] = derivest_of(f, x0, 'DerivativeOrder', 4, 'MethodOrder', 4, 'RombergTerms', 0, 'Style', 'backward')

