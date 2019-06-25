/*

  This file is part of kprime.
  
  kprime is free software: you can redistribute it and/or modify
  it under the terms of the GNU Lesser General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.
  
  kprime is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU Lesser General Public License for more details.
  
  You should have received a copy of the GNU Lesser General Public License
  along with kprime.  If not, see <http://www.gnu.org/licenses/>.


  Created: 2019-06-24
  Copyright: Steven E. Pav, 2019
  Author: Steven E. Pav <shabbychef@gmail.com>
  Comments: Steven E. Pav
*/

#ifndef __DEF_kprime__
#define __DEF_kprime__

#define MAX(_a_,_b_) ((_a_>_b_)? (_a_):(_b_))
#define MIN(_a_,_b_) ((_a_<_b_)? (_a_):(_b_))

#endif /* __DEF_kprime__ */

#include <Rcpp.h>
using namespace Rcpp;

double do_kprime(double q,double v1,double v2,double a,double b,bool lower_tail) {
    // deal with all the corner cases...
    if (a==0) {
        if (b==0) {
            if ((lower_tail && (q >= 0)) || (!lower_tail && (q <= 0))) { return 1.0; } else { return 0.0; }
        } else {
            return Rcpp::pt(q / abs(b),v2,lower_tail,false);
        }
    } else if (b==0) {
        if ((q > 0) && (a < 0)) { return (lower_tail ? 1.0:0.0); }
        if ((q < 0) && (a > 0)) { return (lower_tail ? 0.0:1.0); }
        double anF = (q**2)/(a**2);
        if ((q > 0) && (a > 0)) {
            return Rcpp::pf(anF,v1,v2,lower_tail,false);
        } else {
            return Rcpp::pf(anF,v1,v2,!lower_tail,false);
        }
    } else if (a < 0) {
        // recurse
        return do_kprime(-q,v1,v2,-a,b,!lower_tail);
    } else {
        if (Rcpp::traits::is_infinite<REALSXP>(v1)) {
            if (Rcpp::traits::is_infinite<REALSXP>(v2)) {
                // normal
                return Rcpp::pnorm(q,a,b,lower_tail,false);
            } else {
                // non-central t
                return Rcpp::pnt(q / abs(b),v2,a / abs(b),lower_tail,false);
            }
        } else if (Rcpp::traits::is_infinite<REALSXP>(v2)) {
            // punt on Lambda prime by switching v1 and v2?
            return do_kprime(a,v2,v1,q,b,!lower_tail);
        } else {
            // start here ... 


        }
    }
}

//' @title
//'
//' K-prime distribution
//'
//' @description
//'
//' Distribution, quantile and random variates of the k-prime distribution.
//'
//' @details
//'
//' Suppose \eqn{y \sim \chi^2\left(\nu_1\right)}{y ~ x^2(v1)}, and
//' \eqn{x \sim t \left(\nu_2, a\sqrt{y/\nu_1}/b\right)}{x ~ t(v2,(a/b) sqrt(y/v1))}.
//' Then the random variable
//' \deqn{T = b x}{T = b x}
//' takes a K prime distribution with parameters 
//' \eqn{\nu_1, \nu_2, a, b}{v1, v2, a, b}. In Lecoutre's terminology,
//' \eqn{T \sim K'_{\nu_1, \nu_2}\left(a, b\right)}{T ~ K'_v1,v2(a,b)}
//'
//' Equivalently, we can think of
//' \deqn{T = \frac{b Z + a \sqrt{\chi^2_{\nu_1} / \nu_1}}{\sqrt{\chi^2_{\nu_2} / \nu_2}}}{T = (bZ + a sqrt(chi2_v1/v1)) / sqrt(chi2_v2/v2)}
//' where \eqn{Z} is a standard normal, and the normal and the (central) chi-squares are
//' independent of each other. When \eqn{a=0}{a=0} we recover
//' a central t distribution; 
//' when \eqn{\nu_1=\infty}{v1=inf} we recover a rescaled non-central t distribution;
//' when \eqn{b=0}{b=0}, we get a rescaled square root of a central F
//' distribution; when \eqn{\nu_2=\infty}{v2=inf}, we recover a 
//' Lambda prime distribution.
//'
//' @references
//'
//' Poitevineau, J., and Lecoutre, B. "Implementing Bayesian predictive procedures: The K-prime and K-square distributions." 
//' Computational Statistics & Data Analysis 54, no. 3 (2010): 724-731.
//' \url{https://arxiv.org/abs/1003.4890}
//'
//' @param q vector of quantiles.
//' @param p vector of probabilities.
//' @param n number of observations. 
//' @param v1 the degrees of freedom in the numerator chisquare. When
//' (positive) infinite, we recover a non-central t 
//' distribution with \code{v2} degrees of freedom and non-centrality
//' parameter \code{a}, scaled by \code{b}.
//' This will be recycled against the \code{q,p,n}.
//' @param v2 the degrees of freedom in the denominator chisquare.
//' When equal to infinity, we recover the Lambda prime distribution.
//' This will be recycled against the \code{q,p,n}.
//' @param a the non-centrality scaling parameter. When equal to zero,
//' we recover the (central) t distribution.
//' This will be recycled against the \code{q,p,n}.
//' @param b the scaling parameter.
//' This will be recycled against the \code{q,p,n}.
//' @param lower_tail  boolean of whether to compute the lower tail of the
//' CDF or quantile.
//'
//' @return \code{pkprime} gives the 
//' distribution function, \code{qkprime} gives the quantile function, 
//' and \code{rkprime} generates random deviates.
//'
//' Invalid arguments will result in return value \code{NaN} with a warning.
//'
//' @template etc
//' @name kprime 
//' @rdname kprime 
//' @export
// [[Rcpp::export]]
NumericVector pkprime(NumericVector q, NumericVector v1, NumericVector v2, NumericVector a, NumericVector b, bool lower_tail=true) {
    int lenq = q.length();
    int lenv1 = v1.length();
    int lenv2 = v2.length();
    int lena = a.length();
    int lenb = b.length();
    // sigh
    int maxlen;
    maxlen = MAX(lenq,lenv1);
    maxlen = MAX(maxlen,lenv2);
    maxlen = MAX(maxlen,lena);
    maxlen = MAX(maxlen,lenb);

    NumericVector output(maxlen);
    for (int iii=0;iii < maxlen;iii++) {
        // recycle as needed.
        output[iii] = do_kprime(q[iii % lenq],v1[iii % lenv1],v2[iii % lenv2],a[iii % lena],b[iii % lenb],lower_tail);
    }
    return output;
}


//for vim modeline: (do not edit)
// vim:et:nowrap:ts=4:sw=4:tw=129:fdm=marker:fmr=FOLDUP,UNFOLD:cms=//%s:tags=.c_tags;:syn=cpp:ft=cpp:mps+=<\:>:ai:si:cin:nu:fo=croql:cino=p0t0c5(0:
