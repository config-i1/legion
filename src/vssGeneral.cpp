#include <RcppArmadillo.h>
#include <iostream>
#include <cmath>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

arma::vec vFittedValue(arma::sp_mat const &matrixW, arma::vec const &matrixV, char const &E){
    arma::vec returnedValue;
    switch(E){
        case 'A':
        case 'M':
            returnedValue = matrixW * matrixV;
            break;
        case 'L':
            arma::vec vecYFitted = exp(matrixW * matrixV);
            returnedValue = vecYFitted / sum(vecYFitted);
    }
    return returnedValue;
}

arma::vec vErrorValue(arma::vec const &vectorY, arma::vec const &vectorYFit, char const &E){
    arma::vec returnedValue;
    switch(E){
        case 'A':
        case 'M':
            returnedValue = vectorY - vectorYFit;
        break;
        case 'L':
            arma::vec vectorE = (1 + vectorY - vectorYFit)/2;
            returnedValue = log(vectorE / vectorE(0));
    }
    return returnedValue;
}

/* # Wrapper for fitter */
// [[Rcpp::export]]
arma::cx_vec discounter(arma::sp_mat const &matrixF, arma::sp_mat &matrixW, arma::sp_mat const &matrixG, int const &k){
    arma::cx_vec eigval;
    // If eigen decomposition works, return the values
    // Otherwise return a large number
    if(!arma::eigs_gen(eigval, matrixF - matrixG * matrixW, k)){
        eigval.fill(1e+300);
    };
    return eigval;
    // return arma::conv_to<arma::mat>::from(matrixF - matrixG * matrixW);
}

// Fitter for vector models
List vFitter(arma::mat const &matrixY, arma::mat &matrixV, arma::sp_mat const &matrixF,
             arma::sp_mat &matrixW, arma::sp_mat const &matrixG,
             arma::uvec &lags, char const &E, char const &T, char const &S,
             arma::sp_mat const &matrixO, bool const &backcast,
             unsigned int const &nComponentsTrend) {
    /* matrixY has nrow = nSeries, ncol = obs
     * matrixV has nrow = nSeries * nComponents, ncol = obs + maxlag
     * matrixW, matrixF, matrixG are nSeries * nComponents x nSeries * nComponents.
     * lags is a vector of lags of length nSeries * nComponents
     * matrixX and matrixA are not defined yet.
     */

    unsigned int obs = matrixY.n_cols;
    unsigned int nSeries = matrixY.n_rows;
    unsigned int obsall = matrixV.n_cols;
    // unsigned int nComponents = matrixV.n_rows / nSeries;
    int maxlag = max(lags);
    unsigned int lagsLength = lags.n_rows;
    arma::uvec lagsModifier = lags;
    arma::uvec lagsInternal = lags;

    // Modify lags to select specific cells of matVt
    lagsInternal = lagsInternal * lagsLength;

    for(unsigned int i=0; i<lagsLength; i=i+1){
        lagsModifier(i) = lagsLength - i - 1;
    }

    arma::uvec lagrows(lagsLength, arma::fill::zeros);

    arma::mat matrixYfit(nSeries, obs, arma::fill::zeros);
    arma::mat matrixE(nSeries, obs, arma::fill::zeros);
    // arma::mat bufferforat(matrixGX.n_rows);

    unsigned int const nIterations = 2;

    if(E=='L'){
        matrixW.row(0).zeros();
    }

    // Loop for the backcast
    for (unsigned int j=1; j<=nIterations; j=j+1) {

        /* # Fill in (refine) the head of the matrices */
        for (int i=0; i<maxlag; i=i+1) {
            lagrows = (i+maxlag+1) * lagsLength - (lagsInternal + lagsModifier) - 1;
            matrixV.col(i) = matrixF * matrixV(lagrows);
            // This is needed to reproduce the head properly
            matrixV.col(i+maxlag) = matrixV.col(i);
        }

        for (unsigned int i=maxlag; i<obs+maxlag; i=i+1) {
            lagrows = (i+1) * lagsLength - (lagsInternal + lagsModifier) - 1;

            /* # Measurement equation and the error term */
            matrixYfit.col(i-maxlag) = matrixO.col(i-maxlag) % vFittedValue(matrixW, matrixV(lagrows), E);
            matrixE.col(i-maxlag) = vErrorValue(matrixY.col(i-maxlag), matrixYfit.col(i-maxlag), E);
            // Substitute inf with zero. This might happen for occurrence model
            matrixE.elem(find_nonfinite(matrixE)).fill(0);

            /* # Transition equation */
            matrixV.col(i) = matrixF * matrixV(lagrows) + matrixG * matrixE.col(i-maxlag);
        }

        /* Fill in the tail in the state matrix */
        for (unsigned int i=obs+maxlag; i<obsall; i=i+1) {
            lagrows = (i+1) * lagsLength - (lagsInternal + lagsModifier) - 1;
            matrixV.col(i) = matrixF * matrixV(lagrows);
        }

        ////// Backwards run
        if(backcast && j<(nIterations)){
            // Change the specific element in the state vector to negative
            if(T!='N'){
                matrixV.submat(nSeries,obs+maxlag,nSeries*2-1,obs+maxlag) =
                    -matrixV.submat(nSeries,obs+maxlag,nSeries*2-1,obs+maxlag);
            }

            for (int i=obs+maxlag-1; i>=maxlag; i=i-1) {
                lagrows = (i+1) * lagsLength + lagsInternal - lagsModifier - 1;

                /* # Measurement equation and the error term */
                matrixYfit.col(i-maxlag) = matrixO.col(i-maxlag) % vFittedValue(matrixW, matrixV(lagrows), E);
                matrixE.col(i-maxlag) = vErrorValue(matrixY.col(i-maxlag), matrixYfit.col(i-maxlag), E);
                // Substitute inf with zero. This might happen for occurrence model
                matrixE.elem(find_nonfinite(matrixE)).fill(0);

                /* # Transition equation */
                matrixV.col(i) = matrixF * matrixV(lagrows) + matrixG * matrixE.col(i-maxlag);
            }

            /* # Fill in the head of the matrices */
            for (int i=maxlag-1; i>=0; i=i-1) {
                lagrows = (i+1) * lagsLength + lagsInternal - lagsModifier - 1;
                matrixV.col(i) = matrixF * matrixV(lagrows);
            }

            // Change the specific element in the state vector to negative
            if(T!='N'){
                matrixV.submat(nSeries,0,nSeries*2-1,0) =
                    -matrixV.submat(nSeries,0,nSeries*2-1,0);
            }
        }
    }

    // , Named("matat") = matrixA
    return List::create(Named("matVt") = matrixV, Named("yfit") = matrixYfit,
                        Named("errors") = matrixE);
}

/* # Wrapper for fitter */
// [[Rcpp::export]]
RcppExport SEXP vFitterWrap(arma::mat const &matrixY, arma::mat matrixV, arma::sp_mat &matrixF,
                            arma::sp_mat &matrixW, arma::sp_mat &matrixG,
                            arma::uvec &lags, char const &E, char const &T, char const &S,
                            arma::sp_mat &matrixO, bool const &backcast,
                            unsigned int const &nComponentsTrend) {
// SEXP matxt, SEXP matat, SEXP matFX, SEXP matGX,
    // NumericMatrix yt_n(yt);
    // arma::mat matrixY(yt_n.begin(), yt_n.nrow(), yt_n.ncol(), false);

    // NumericMatrix matVt_n(matVt);
    // arma::mat matrixV(matVt_n.begin(), matVt_n.nrow(), matVt_n.ncol());

    // NumericMatrix matF_n(matF);
    // arma::mat matrixF(matF_n.begin(), matF_n.nrow(), matF_n.ncol(), false);

    // NumericMatrix matw_n(matw);
    // arma::mat matrixW(matw_n.begin(), matw_n.nrow(), matw_n.ncol(), false);

    // NumericMatrix matG_n(matG);
    // arma::mat matrixG(matG_n.begin(), matG_n.nrow(), matG_n.ncol(), false);

    // IntegerVector modellags_n(modellags);
    // arma::uvec lags = as<arma::uvec>(modellags_n);

    // char E = as<char>(Etype);
    // char T = as<char>(Ttype);
    // char S = as<char>(Stype);

    // NumericMatrix matxt_n(matxt);
    // arma::mat matrixX(matxt_n.begin(), matxt_n.nrow(), matxt_n.ncol(), false);
    //
    // NumericMatrix matat_n(matat);
    // arma::mat matrixA(matat_n.begin(), matat_n.nrow(), matat_n.ncol());
    //
    // NumericMatrix matFX_n(matFX);
    // arma::mat matrixFX(matFX_n.begin(), matFX_n.nrow(), matFX_n.ncol(), false);
    //
    // NumericMatrix matGX_n(matGX);
    // arma::mat matrixGX(matGX_n.begin(), matGX_n.nrow(), matGX_n.ncol(), false);

    // NumericMatrix ot_n(ot);
    // arma::mat matrixO(ot_n.begin(), ot_n.nrow(), ot_n.ncol(), false);

    return wrap(vFitter(matrixY, matrixV, matrixF,
                        matrixW, matrixG,
                        lags, E, T, S,
                        matrixO, backcast,
                        nComponentsTrend));
}


/* # Function produces the point forecasts for the specified model */
arma::mat vForecaster(arma::mat const & matrixV, arma::sp_mat const &matrixF, arma::sp_mat matrixW,
                      unsigned int const &nSeries, unsigned int const &hor,
                      char const &E, char const &T, char const &S, arma::uvec lags){
                      // arma::mat const &matrixX, arma::mat const &matrixA, arma::mat const &matrixFX
    int lagsLength = lags.n_rows;
    unsigned int maxlag = max(lags);
    unsigned int hh = hor + maxlag;

    arma::uvec lagrows(lagsLength, arma::fill::zeros);
    arma::mat matYfor(nSeries, hor, arma::fill::zeros);
    arma::mat matrixVnew(matrixV.n_rows, hh, arma::fill::zeros);
    // arma::mat matrixAnew(hh, matrixA.n_cols, arma::fill::zeros);

    lags = lags * lagsLength;

    for(int i=0; i<lagsLength; i=i+1){
        lags(i) = lags(i) + (lagsLength - i - 1);
    }

    if(E=='L'){
        matrixW.row(0).zeros();
    }

    matrixVnew.submat(0,0,matrixVnew.n_rows-1,maxlag-1) = matrixV.submat(0,0,matrixVnew.n_rows-1,maxlag-1);
    // matrixAnew.submat(0,0,maxlag-1,matrixAnew.n_cols-1) = matrixAnew.submat(0,0,maxlag-1,matrixAnew.n_cols-1);

    /* # Fill in the new xt matrix using F. Do the forecasts. */
    for (unsigned int i=maxlag; i<hh; i=i+1) {
        lagrows = (i+1) * lagsLength - lags - 1;

        /* # Transition equation */
        matrixVnew.col(i) = matrixF * matrixVnew(lagrows);
        // matrixAnew.row(i) = matrixAnew.row(i-1) * matrixFX;

        matYfor.col(i-maxlag) = vFittedValue(matrixW, matrixVnew(lagrows), E);
        // matYfor.col(i-maxlag) = matrixW * matrixVnew(lagrows);
    }

    return matYfor;
}

/* # Wrapper for forecaster */
// [[Rcpp::export]]
RcppExport SEXP vForecasterWrap(arma::mat matrixV, arma::sp_mat const &matrixF, arma::sp_mat const &matrixW,
                                unsigned int const &nSeries, unsigned int const &hor,
                                char const &E, char const &T, char const &S, arma::uvec &lags){
    // SEXP matxt, SEXP matat, SEXP matFX

    // NumericMatrix matVt_n(matVt);
    // arma::mat matrixV(matVt_n.begin(), matVt_n.nrow(), matVt_n.ncol(), false);

    // NumericMatrix matF_n(matF);
    // arma::mat matrixF(matF_n.begin(), matF_n.nrow(), matF_n.ncol(), false);
    //
    // NumericMatrix matw_n(matw);
    // arma::mat matrixW(matw_n.begin(), matw_n.nrow(), matw_n.ncol(), false);

    // unsigned int nSeries = as<int>(series);
    // unsigned int hor = as<int>(h);
    // char E = as<char>(Etype);
    // char T = as<char>(Ttype);
    // char S = as<char>(Stype);
    //
    // IntegerVector modellags_n(modellags);
    // arma::uvec lags = as<arma::uvec>(modellags_n);

    // NumericMatrix matxt_n(matxt);
    // arma::mat matrixX(matxt_n.begin(), matxt_n.nrow(), matxt_n.ncol(), false);
    //
    // NumericMatrix matat_n(matat);
    // arma::mat matrixA(matat_n.begin(), matat_n.nrow(), matat_n.ncol());
    //
    // NumericMatrix matFX_n(matFX);
    // arma::mat matrixFX(matFX_n.begin(), matFX_n.nrow(), matFX_n.ncol(), false);

    return wrap(vForecaster(matrixV, matrixF, matrixW, nSeries, hor, E, T, S, lags));
}

/* # Function returns the chosen Cost Function based on the chosen model and produced errors */
// double vOptimiser(arma::mat const &matrixY, arma::mat &matrixV, arma::sp_mat const &matrixF,
//                   arma::sp_mat &matrixW, arma::sp_mat const &matrixG,
//                   arma::uvec &lags, char const &E, char const &T, char const &S,
//                   char const& CFtype, double const &normalize, arma::sp_mat const &matrixO, arma::mat matrixOtObs){
//     // bool const &multi, std::string const &CFtype, char const &fitterType,
//     // arma::mat const &matrixX, arma::mat &matrixA, arma::mat const &matrixFX, arma::mat const &matrixGX,
//     // # Make decomposition functions shut up!
//     // std::ostream nullstream(0);
//     // arma::arma_cerr_stream<char>(&nullstream);
//
//     arma::uvec nonzeroes = find(matrixO>0);
//     int obs = nonzeroes.n_rows;
//     double CFres = 0;
//
//     int nSeries = matrixY.n_rows;
//
//     List fitting = vFitter(matrixY, matrixV, matrixF, matrixW, matrixG, lags, E, T, S, matrixO);
//
//     NumericMatrix mvtfromfit = as<NumericMatrix>(fitting["matVt"]);
//     matrixV = as<arma::mat>(mvtfromfit);
//     NumericMatrix errorsfromfit = as<NumericMatrix>(fitting["errors"]);
//
//     arma::mat matErrors(errorsfromfit.begin(), errorsfromfit.nrow(), errorsfromfit.ncol(), false);
//     // matErrors = matErrors / matrixOtObs;
//
//     if(E=='L'){
//         NumericMatrix Yfromfit = as<NumericMatrix>(fitting["yfit"]);
//         arma::mat matrixYfit(Yfromfit.begin(), Yfromfit.nrow(), Yfromfit.ncol(), false);
//         CFres = -sum(log(matrixYfit.elem(arma::find(matrixY==1))));
//     }
//     else{
//         if(CFtype=='l'){
//             try{
//                 CFres = double(log(arma::prod(eig_sym((matErrors / normalize) * arma::trans(matErrors / normalize) / matrixOtObs))) +
//                     nSeries * log(pow(normalize,2)));
//             }
//             catch(const std::runtime_error&){
//                 CFres = double(log(arma::det((matErrors / normalize) * arma::trans(matErrors / normalize) / matrixOtObs)) +
//                     nSeries * log(pow(normalize,2)));
//             }
//         }
//         else if(CFtype=='d'){
//             CFres = arma::as_scalar(sum(log(sum(pow(matErrors,2)) / double(obs)), 1));
//         }
//         else{
//             CFres = arma::as_scalar(sum(sum(pow(matErrors,2)) / double(obs), 1));
//         }
//     }
//     return CFres;
// }
//
//
// /* # This is a wrapper for optimizer, which currently uses admissible bounds */
// // [[Rcpp::export]]
// RcppExport SEXP vOptimiserWrap(arma::mat const &matrixY, arma::mat matrixV,
//                                arma::sp_mat const &matrixF, arma::sp_mat &matrixW, arma::sp_mat &matrixG,
//                                arma::uvec &lags, char const &E, char const &T, char const &S,
//                                char const &CFtype, double const &normalize, char const &boundtype,
//                                arma::sp_mat const &matrixO, arma::mat const &matrixOtObs) {
//     // SEXP multisteps, SEXP CFt, SEXP fittertype, SEXP bounds,
//     // SEXP matxt, SEXP matat, SEXP matFX, SEXP matGX
//     /* Function is needed to implement admissible constrains on smoothing parameters */
//     // NumericMatrix yt_n(yt);
//     // arma::mat matrixY(yt_n.begin(), yt_n.nrow(), yt_n.ncol(), false);
//
//     // NumericMatrix matVt_n(matVt);
//     // arma::mat matrixV(matVt_n.begin(), matVt_n.nrow(), matVt_n.ncol());
//
//     // NumericMatrix matF_n(matF);
//     // arma::mat matrixF(matF_n.begin(), matF_n.nrow(), matF_n.ncol(), false);
//
//     // NumericMatrix matw_n(matw);
//     // arma::mat matrixW(matw_n.begin(), matw_n.nrow(), matw_n.ncol(), false);
//     //
//     // NumericMatrix matG_n(matG);
//     // arma::mat matrixG(matG_n.begin(), matG_n.nrow(), matG_n.ncol(), false);
//     //
//     // IntegerVector modellags_n(modellags);
//     // arma::uvec lags = as<arma::uvec>(modellags_n);
//
//     // char E = as<char>(Etype);
//     // char T = as<char>(Ttype);
//     // char S = as<char>(Stype);
//
//     // char CFtype = as<char>(cfType);
//
//     // char fitterType = as<char>(fittertype);
//
//     // char boundtype = as<char>(bounds);
//
//     // double normalize = as<double>(normalizer);
//
//     // NumericMatrix matxt_n(matxt);
//     // arma::mat matrixX(matxt_n.begin(), matxt_n.nrow(), matxt_n.ncol(), false);
//     //
//     // NumericMatrix matat_n(matat);
//     // arma::mat matrixA(matat_n.begin(), matat_n.nrow(), matat_n.ncol());
//     //
//     // NumericMatrix matFX_n(matFX);
//     // arma::mat matrixFX(matFX_n.begin(), matFX_n.nrow(), matFX_n.ncol(), false);
//     //
//     // NumericMatrix matGX_n(matGX);
//     // arma::mat matrixGX(matGX_n.begin(), matGX_n.nrow(), matGX_n.ncol(), false);
//
//     // NumericMatrix ot_n(ot);
//     // arma::mat matrixO(ot_n.begin(), ot_n.nrow(), ot_n.ncol(), false);
//
//     // NumericMatrix otObs_n(otObs);
//     // arma::mat matrixOtObs(otObs_n.begin(), otObs_n.nrow(), otObs_n.ncol(), false);
//
//     // Values needed for eigenvalues calculation
//     arma::cx_vec eigval;
//
//     if(boundtype=='a'){
//         if(arma::eig_gen(eigval, matrixF - matrixG * matrixW)){
//             if(max(abs(eigval)) > (1 + 1E-50)){
//                 return wrap(max(abs(eigval))*1E+100);
//             }
//         }
//         else{
//             return wrap(1E+300);
//         }
//     }
//
//     // multi, CFtype, fitterType,
//     return wrap(vOptimiser(matrixY, matrixV, matrixF, matrixW, matrixG,
//                            lags, E, T, S, CFtype, normalize, matrixO, matrixOtObs));
// }
