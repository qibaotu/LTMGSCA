#include<cstdio>
#include<cstdlib>
#include<cmath>
#include<iostream>
#include<algorithm>
#include<fstream>
#include<Rcpp.h>
#include<list>

double cdf(double x, double m, double s)
{
  return 0.5 * (1 + erf((x - m) / (s * sqrt(2.))));
}

double pdf(double x, double m, double s) {
  static const double inv_sqrt_2pi = 0.3989422804014327;
  double a = (x - m) / s;
  return inv_sqrt_2pi / s * std::exp(-0.5f * a * a);
}

Rcpp::NumericVector Hf(const double Zcut, const Rcpp::NumericVector ul, const Rcpp::NumericVector sigmal, const int k) {
  double a, b, c;
  Rcpp::NumericVector d;
  for (int i = 0; i < k; i++) {
    c = (Zcut - ul[i]) / sigmal[i];
    a = pdf(c, 0, 1);
    b = cdf(c, 0, 1);
    if (b != 0) {
      d.push_back(a / b);
    }
    else {
      d.push_back(a);
    }
  }
  return d;
}

Rcpp::NumericVector Pi_Zj_Zcut_new(const double Zcut, const Rcpp::NumericVector ul, const Rcpp::NumericVector sigmal, const Rcpp::NumericVector wl0, const int k) {
  Rcpp::NumericVector a(k);
  double sum = 0.0;
  for (int i = 0; i < k; i++) {
    a[i] = wl0[i] * cdf(Zcut, ul[i], sigmal[i]);
    sum += a[i];
  }
  for (int i = 0; i < k; i++) {
    a[i] = a[i] / sum;
  }
  return a;
}

double sigma_f(const Rcpp::NumericVector wl, const Rcpp::NumericVector ul, const Rcpp::NumericVector sigmal, const Rcpp::NumericVector y, const int i, const int j, const int k) {
  double sigma_t = 0.0;
  sigma_t = wl[j] * pdf(y[i], ul[j], sigmal[j]);
  double temp = 0.0;
  for (int id = 0; id < k; id++) {
    temp += wl[id] * pdf(y[i], ul[id], sigmal[id]);
  }
  sigma_t /= temp;
  return sigma_t;
}
// [[Rcpp::export]]
Rcpp::List Separate_K_rpkm_new1(Rcpp::NumericVector y, const int ROUNDS, const int Zcut, const int k, double err = 1e-10) {
  double Cutted = 0.0;
  Rcpp::NumericVector ybig;
  int ybig_id_num = 0;
  double ybig_sum = 0.0;
  for (auto it = y.begin(); it != y.end(); it++) {
    if (*it < Zcut) {
      Cutted++;
    }
    else {
      ybig.push_back(*it);
      ybig_id_num++;
      ybig_sum += *it;
    }
  }
  double ybig_mean = ybig_sum / ybig_id_num;
  double temp_sum = 0.0;
  for (int i = 0; i < ybig_id_num; i++) {
    temp_sum += (ybig[i] - ybig_mean)*(ybig[i] - ybig_mean);
  }
  double var = temp_sum / (ybig_id_num - 1);
  double sig = sqrt(var);
  Rcpp::NumericVector wl1(k);
  Rcpp::NumericVector wl(k);
  Rcpp::NumericVector ul1(k);
  Rcpp::NumericVector ul(k);
  Rcpp::NumericVector sigmal1(k);
  Rcpp::NumericVector sigmal(k);
  int n = ybig_id_num;
  Rcpp::NumericVector ul0;
  std::sort(ybig.begin(), ybig.end());
  for (int i = 1; i <= k; i++) {
    int temp = floor(i*n / (k + 1));
    ul0.push_back(ybig[temp - 1]);
  }
  ul0[0] = *std::min_element(ybig.begin(), ybig.end()) - 1;
  ul0[k - 1] = *std::max_element(ybig.begin(), ybig.end()) + 1;
  Rcpp::NumericVector wl0(k);
  Rcpp::NumericVector sigmal0(k);
  for (int i = 0; i < k; i++) {
    wl0[i] = (double)1 / k;
    sigmal0[i] = sig;
  }
  for (int i = 0; i < k; i++) {
    wl[i] = wl0[i];				//wl<-wl0
    ul[i] = ul0[i];				//ul<-ul0
    sigmal[i] = sigmal0[i];		//sigmal<-sigmal0
  }
  Rcpp::NumericMatrix table_theta_t1(ybig_id_num, k);
  /*	for (int i = 0; i < table_theta_t1.nrow(); i++)
  for (int j = 0; j < table_theta_t1.ncol(); j++)
  table_theta_t1(i,j) = 0.0;
  */
  for (int rounds = 0; rounds < ROUNDS; rounds++) {
    for (int i = 0; i < table_theta_t1.nrow(); i++) {
      for (int j = 0; j < table_theta_t1.ncol(); j++) {
        table_theta_t1(i,j) = sigma_f(wl, ul, sigmal, ybig, i, j, k);
      }
    }
    Rcpp::NumericVector pZil;
    pZil = Pi_Zj_Zcut_new(Zcut, ul, sigmal, wl, k);
    Rcpp::NumericVector table_theta_t1_sum;
    Rcpp::NumericVector denoml;
    double wl1_sum = 0.0;
    for (int i = 0; i < k; i++) {
      table_theta_t1_sum.push_back(0.0);
      for (int j = 0; j < ybig_id_num; j++) {
        table_theta_t1_sum[i] += table_theta_t1(j,i);
      }
      denoml.push_back(table_theta_t1_sum[i] + pZil[i] * Cutted);
      wl1[i] = denoml[i] / (ybig_id_num + Cutted);
      wl1_sum += wl1[i];
    }
    for (int i = 0; i < k; i++) {
      wl1[i] = wl1[i] / wl1_sum;
    }
    Rcpp::NumericVector hf;
    hf = Hf(Zcut, ul, sigmal, k);
    for (int id = 0; id < k; id++) {
      double ybig_t1_sum = 0.0;
      for (int i = 0; i < ybig_id_num; i++) {
        ybig_t1_sum += ybig[i] * table_theta_t1(i,id);
      }
      ul1[id] = (ybig_t1_sum + (ul[id] - sigmal[id] * hf[id])*pZil[id] * Cutted) / denoml[id];
    }
    for (int id = 0; id < k; id++) {
      double ybig_ul_t1_sum = 0.0;
      for (int i = 0; i < ybig_id_num; i++) {
        ybig_ul_t1_sum += (ybig[i] - ul[id]) * (ybig[i] - ul[id]) * table_theta_t1(i,id);
      }

      sigmal1[id] = sqrt((ybig_ul_t1_sum + sigmal[id] * sigmal[id]
                            * (1 - (Zcut - ul[id]) / sigmal[id] * hf[id]) * pZil[id] * Cutted) / denoml[id]);
    }
    double wlabs = 0;
    double ulabs = 0;
    double sigmalabs = 0;
    for(int i = 0; i < k; i++){
      wlabs += fabs(wl[i] - wl1[i]);
      ulabs += fabs(ul[i] - ul1[i]);
      sigmalabs += fabs(sigmal[i] - sigmal1[i]);
    }
    wlabs /= k;
    ulabs /= k;
    sigmalabs /= k;
    if(wlabs <= err && ulabs <= err && sigmalabs <= err){
      std::cout << "rounds = " << rounds << std::endl;
      break;
    }

    wl.assign(wl1.begin(), wl1.end());
    ul.assign(ul1.begin(), ul1.end());
    sigmal.assign(sigmal1.begin(), sigmal1.end());
  }
  return Rcpp::List::create(
    Rcpp::Named("wlValue") = wl,
    Rcpp::Named("ulValue") = ul,
    Rcpp::Named("sigmalValue") = sigmal);
}
