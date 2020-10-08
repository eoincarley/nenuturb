function fit_psd, frequency, power, pspecerr=pspecerr

        start = [-1, -1]
        fit = 'p[0] + p[1]*x'
        err = power
        err[*] = pspecerr*abs(power)

        ; Note here that if I supply the weights keyword with poisson weighting the fit
        ; gives a p-value of 1 for eveything, e.g. perfect fit within error. Every fit
        ; is then accepted below, even a few power spectra that are not power law.
        ; If I give an errors (without the weights) then not everything is accepted. Some
        ; fits are rejected. However, the choice of err here is slightly arbitrary.

        p = mpfitexpr(fit, frequency, power, err, weights=1/power^2, start, perror = perror, $
                yfit=yfit, bestnorm=bestnorm, dof=dof, /quiet)
        perror = perror * SQRT(BESTNORM / DOF)
        aerr = perror[1]*2.0 ; 2-sigma uncertainty on the slope
        ierr = perror[0]*2.0 ; 2-sigma uncertainty on the intercept

        ;----------------------------------------;
        ; This is the probability the chi^2 score 
        ; is worse than calculated value.
        ; If <5% reject the fit
        pvalue = (1.0-CHISQR_PDF(bestnorm, DOF))*100.0

        return, [p[0], p[1], aerr, ierr, pvalue]

end
