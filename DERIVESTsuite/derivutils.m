classdef derivutils
    methods(Static)
    function make_fdarules()
        StepRatio = 2.0000001;
        Styles = {'Central', 'Forward', 'Backward'};
        fprintf('double* get_fdarule(int derivative_order, int method_order, DerivestStyle style, int *n) {\n');
        for style_index = 1:3
            Style = Styles{style_index};
            for MethodOrder=1:4
                for DerivativeOrder=1:4
                    if (lower(Style(1))=='c') && (mod(MethodOrder,2)==1), continue; end
                    fdarule = derivutils.get_fdarule(Style, MethodOrder, DerivativeOrder, StepRatio);
                    fprintf('  if (style==DerivestStyle_%s && method_order==%i && derivative_order==%i) { static double data[] = { %s}; *n = %i; return data; }\n', Style, MethodOrder, DerivativeOrder, sprintf('%.20f, ', fdarule), length(fdarule));
                end
            end
        end
        fprintf('  return NULL;\n');
        fprintf('}\n');
    end

    function make_rombmat()
        StepRatio = 2.0000001;
        func = {'double* get_qromb(DerivestStyle style, int method_order, int romberg_terms, double* err, int *rows, int *cols) {\n', ...
                'double* get_rmat(DerivestStyle style, int method_order, int romberg_terms, double* err, int *rows, int *cols) {\n', ...
                'double* get_rinv(DerivestStyle style, int method_order, int romberg_terms, double* err, int *rows, int *cols) {\n'};

        for func_index=1:3
            fprintf(func{func_index});
            Styles = {'Central', 'Forward', 'Backward'};
            for style_index = 1:3
                Style = Styles{style_index};
                for MethodOrder=1:4
                    for RombergTerms=0:3
                        if (lower(Style(1))=='c') && (mod(MethodOrder,2)==1), continue; end
                        [qromb,rmat,rinv,err] = derivutils.rombmat(StepRatio,lower(Style(1))=='c', RombergTerms, MethodOrder);
                        matrices = {qromb, rmat, rinv}; matrix = matrices{func_index}; matrix_transpose = matrix';%
                        matrix_as_str = sprintf('%.20f, ', matrix_transpose(:));
                        % alternatively, pass values as binary - and don't forget to use "unsigned long long" and "reinterpret_cast<double*>(data)" in C code
                        % matrix_as_str = sprintf('%u, ', typecast(matrix_transpose(:), 'uint64'));
                        fprintf('  if (style==DerivestStyle_%s && method_order==%i && romberg_terms==%i) { static double data[] = { %s}; *rows = %i; *cols = %i; *err=%.20f; return data; }\n', Style, MethodOrder, RombergTerms, matrix_as_str, size(matrix, 1), size(matrix, 2), err);
                    end
                end
            end
            fprintf('  return NULL;\n');
            fprintf('}\n');
        end
    end
    
    function make_c_code()
        derivutils.make_fdarules();
        derivutils.make_rombmat();
    end
    
    function fdarule = get_fdarule(Style, MethodOrder, DerivativeOrder, StepRatio)
        % generate finite differencing rule in advance.
        % The rule is for a nominal unit step size, and will
        % be scaled later to reflect the local step size.
        fdarule = 1;
        switch lower(Style)
          case 'central'
            % for central rules, we will reduce the load by an
            % even or odd transformation as appropriate.
            if MethodOrder==2
              switch DerivativeOrder
                case 1
                  % the odd transformation did all the work
                  fdarule = 1;
                case 2
                  % the even transformation did all the work
                  fdarule = 2;
                case 3
                  % the odd transformation did most of the work, but
                  % we need to kill off the linear term
                  fdarule = [0 1]/derivutils.fdamat(StepRatio,1,2);
                case 4
                  % the even transformation did most of the work, but
                  % we need to kill off the quadratic term
                  fdarule = [0 1]/derivutils.fdamat(StepRatio,2,2);
              end
            else
              % a 4th order method. We've already ruled out the 1st
              % order methods since these are central rules.
              switch DerivativeOrder
                case 1
                  % the odd transformation did most of the work, but
                  % we need to kill off the cubic term
                  fdarule = [1 0]/derivutils.fdamat(StepRatio,1,2);
                case 2
                  % the even transformation did most of the work, but
                  % we need to kill off the quartic term
                  fdarule = [1 0]/derivutils.fdamat(StepRatio,2,2);
                case 3
                  % the odd transformation did much of the work, but
                  % we need to kill off the linear & quintic terms
                  fdarule = [0 1 0]/derivutils.fdamat(StepRatio,1,3);
                case 4
                  % the even transformation did much of the work, but
                  % we need to kill off the quadratic and 6th order terms
                  fdarule = [0 1 0]/derivutils.fdamat(StepRatio,2,3);
              end
            end
          case {'forward' 'backward'}
            % These two cases are identical, except at the very end,
            % where a sign will be introduced.

            % No odd/even trans, but we already dropped
            % off the constant term
            if MethodOrder==1
              if DerivativeOrder==1
                % an easy one
                fdarule = 1;
              else
                % 2:4
                v = zeros(1,DerivativeOrder);
                v(DerivativeOrder) = 1;
                fdarule = v/derivutils.fdamat(StepRatio,0,DerivativeOrder);
              end
            else
              % MethodOrder methods drop off the lower order terms,
              % plus terms directly above DerivativeOrder
              v = zeros(1,DerivativeOrder + MethodOrder - 1);
              v(DerivativeOrder) = 1;
              fdarule = v/derivutils.fdamat(StepRatio,0,DerivativeOrder+MethodOrder-1);
            end

            % correct sign for the 'backward' rule
            if lower(Style(1)) == 'b'
              fdarule = -fdarule;
            end

        end % switch on style (generating fdarule)
    end % get_fdarule
    
    function mat = fdamat(sr,parity,nterms)
        % Compute matrix for fda derivation.
        % parity can be
        %   0 (one sided, all terms included but zeroth order)
        %   1 (only odd terms included)
        %   2 (only even terms included)
        % nterms - number of terms

        % sr is the ratio between successive steps
        srinv = 1./sr;

        switch parity
          case 0
            % single sided rule
            [i,j] = ndgrid(1:nterms);
            c = 1./factorial(1:nterms);
            mat = c(j).*srinv.^((i-1).*j);
          case 1
            % odd order derivative
            [i,j] = ndgrid(1:nterms);
            c = 1./factorial(1:2:(2*nterms));
            mat = c(j).*srinv.^((i-1).*(2*j-1));
          case 2
            % even order derivative
            [i,j] = ndgrid(1:nterms);
            c = 1./factorial(2:2:(2*nterms));
            mat = c(j).*srinv.^((i-1).*(2*j));
        end
      end % fdamat

      function [qromb,rmat,rinv,err] = rombmat(StepRatio,isStyleCentral,RombergTerms,MethodOrder)
        % calculate romberg matrices
        % errest - error estimates
        if isStyleCentral
          rombexpon = 2*(1:RombergTerms) + MethodOrder - 2;
        else
          rombexpon = (1:RombergTerms) + MethodOrder - 1;
        end
        
        srinv = 1/StepRatio;

        % do nothing if no romberg terms
        nexpon = length(rombexpon);
        rmat = ones(nexpon+2,nexpon+1);
        switch nexpon
          case 0
            % rmat is simple: ones(2,1)
          case 1
            % only one romberg term
            rmat(2,2) = srinv^rombexpon;
            rmat(3,2) = srinv^(2*rombexpon);
          case 2
            % two romberg terms
            rmat(2,2:3) = srinv.^rombexpon;
            rmat(3,2:3) = srinv.^(2*rombexpon);
            rmat(4,2:3) = srinv.^(3*rombexpon);
          case 3
            % three romberg terms
            rmat(2,2:4) = srinv.^rombexpon;
            rmat(3,2:4) = srinv.^(2*rombexpon);
            rmat(4,2:4) = srinv.^(3*rombexpon);
            rmat(5,2:4) = srinv.^(4*rombexpon);
        end

        % qr factorization used for the extrapolation as well
        % as the uncertainty estimates
        [qromb,rromb] = qr(rmat,0);

        % the noise amplification is further amplified by the Romberg step.
        % amp = cond(rromb);

        rinv = rromb\eye(nexpon+1);

        % uncertainty estimate of derivative prediction
        cov1 = sum(rinv.^2,2); % 1 spare dof
        err = 12.7062047361747*sqrt(cov1(1));
      end % rombmat
        
      function [der_romb,errest] = rombextrap(der_init, qromb, rmat, rinv, err)
        % do romberg extrapolation for each estimate
        %
        %  StepRatio - Ratio decrease in step
        %  der_init - initial derivative estimates
        %  rombexpon - higher order terms to cancel using the romberg step
        %
        %  der_romb - derivative estimates returned
        %  errest - eror estimates
        %  amp - noise amplification factor due to the romberg step

        % this does the extrapolation to a zero step size.
        ne = length(der_init);
        nexpon = size(rmat, 1) - 2;
        nrombcoefs = max(1,ne - (nexpon+2));
        
        % implementing matrix product as a custom code gives slightly different result (due to numerical issues)
        % set this flag to 1 to reproduce exact DERIVESTsuite results.
        % set this flag to 0 to reproduce cppDERIVEST results.
        % the differences between these two are tiny
        use_DERIVESTsuite_numerics = 1;

        if use_DERIVESTsuite_numerics
            rhs = derivutils.vec2mat(der_init,nexpon+2,max(1,ne - (nexpon+2)));
            rombcoefs = rinv * (qromb.'*rhs); 
        else 
            rombcoefs = zeros(size(rinv, 2), nrombcoefs);
            for i=1:size(rinv, 1)
                for j=1:nrombcoefs
                    for k=1:size(rinv, 2)
                        for q = 1:size(qromb, 1)
                            rombcoefs(i,j)=rombcoefs(i,j)+rinv(i,k) * qromb(q, k) * der_init(q + j - 1);
                        end
                        %another alternative is to first multiply qromb and rhs
                        %qromb_rhs = zeros(size(qromb, 2), size(rhs, 2));
                        %for k=1:size(qromb, 2), for j=1:size(rhs, 2), for q=1:size(qromb, 1)
                        %    qromb_rhs(k, j) = qromb_rhs(k,j) + qromb(q, k)*der_init(q + j - 1);  % qromb(k, i)*rhs(k, j);
                        %end; end; end                    
                        %rombcoefs(i,j)=rombcoefs(i,j)+rinv(i,k) * qromb_rhs(k, j);
                    end
                end
            end
        end
        
        der_romb = rombcoefs(1,:).';

        % uncertainty estimate of derivative prediction
        if use_DERIVESTsuite_numerics
            s = sqrt(sum((rhs - rmat*rombcoefs).^2,1));
        else
            s = zeros(1, nrombcoefs);
            for j=1:nrombcoefs
                for i=1:size(rmat, 1)
                    rmat_rombcoefs = 0;
                    for k=1:size(rmat, 2)
                        rmat_rombcoefs=rmat_rombcoefs+rmat(i,k)*rombcoefs(k,j); 
                    end
                    s(j) = s(j) + (der_init(i + j - 1) - rmat_rombcoefs).^2;
                end
                s(j) = sqrt(s(j));
            end
        end
        errest = s.'*err;

      end % rombextrap
      function mat = vec2mat(vec,n,m)
        % forms the matrix M, such that M(i,j) = vec(i+j-1)
        [i,j] = ndgrid(1:n,0:m-1);
        ind = i+j;
        mat = vec(ind);
        if n==1
          mat = mat.';
        end
      end % vec2mat
   end % methods(Static)
end % classdef derivutils
