// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#ifndef RCPP_bvartools_RCPPEXPORTS_H_GEN_
#define RCPP_bvartools_RCPPEXPORTS_H_GEN_

#include <RcppArmadillo.h>
#include <Rcpp.h>

namespace bvartools {

    using namespace Rcpp;

    namespace {
        void validateSignature(const char* sig) {
            Rcpp::Function require = Rcpp::Environment::base_env()["require"];
            require("bvartools", Rcpp::Named("quietly") = true);
            typedef int(*Ptr_validate)(const char*);
            static Ptr_validate p_validate = (Ptr_validate)
                R_GetCCallable("bvartools", "_bvartools_RcppExport_validate");
            if (!p_validate(sig)) {
                throw Rcpp::function_not_exported(
                    "C++ function with signature '" + std::string(sig) + "' not found in bvartools");
            }
        }
    }

    inline arma::vec stoch_vol(arma::vec y, arma::vec h, double sigma, double h_init) {
        typedef SEXP(*Ptr_stoch_vol)(SEXP,SEXP,SEXP,SEXP);
        static Ptr_stoch_vol p_stoch_vol = NULL;
        if (p_stoch_vol == NULL) {
            validateSignature("arma::vec(*stoch_vol)(arma::vec,arma::vec,double,double)");
            p_stoch_vol = (Ptr_stoch_vol)R_GetCCallable("bvartools", "_bvartools_stoch_vol");
        }
        RObject rcpp_result_gen;
        {
            RNGScope RCPP_rngScope_gen;
            rcpp_result_gen = p_stoch_vol(Shield<SEXP>(Rcpp::wrap(y)), Shield<SEXP>(Rcpp::wrap(h)), Shield<SEXP>(Rcpp::wrap(sigma)), Shield<SEXP>(Rcpp::wrap(h_init)));
        }
        if (rcpp_result_gen.inherits("interrupted-error"))
            throw Rcpp::internal::InterruptedException();
        if (Rcpp::internal::isLongjumpSentinel(rcpp_result_gen))
            throw Rcpp::LongjumpException(rcpp_result_gen);
        if (rcpp_result_gen.inherits("try-error"))
            throw Rcpp::exception(Rcpp::as<std::string>(rcpp_result_gen).c_str());
        return Rcpp::as<arma::vec >(rcpp_result_gen);
    }

}

#endif // RCPP_bvartools_RCPPEXPORTS_H_GEN_
