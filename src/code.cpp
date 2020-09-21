#include <Rcpp.h>
#include <RcppEigen.h>
using namespace Rcpp;
// using namespace Eigen;


// [[Rcpp::depends(RcppEigen)]]

// // [[Rcpp::export]]
// NumericMatrix cov_cpp(NumericMatrix x,
//                       NumericMatrix y){
//   int n_x = x.cols();
//   int n_y = y.cols();
//   int n_rows = x.rows();
//
//   NumericMatrix out_mat(n_x, n_y);
//
//   for(int i=0; i<n_x; i++){
//     for(int j=0; j<n_y; j++){
//       double mu_x = sum(x.column(i)) / n_rows;
//       double mu_y = sum(y.column(j)) / n_rows;
//       double s=0;
//
//       for(int k=0; k<n_rows; k++){
//         s += (x(k, i) - mu_x) * (y(k, j) - mu_y);
//       }
//
//       out_mat(i, j) = s / (n_rows - 1);
//     }
//   }
//   return(out_mat);
// }

// [[Rcpp::export]]
Eigen::MatrixXd cov_cpp2(Eigen::MatrixXd x,
                         Eigen::MatrixXd y){
  int n_x = x.cols();
  int n_y = y.cols();
  int n_rows = x.rows();

  Eigen::MatrixXd out_mat(n_x, n_y);

  for(int i=0; i<n_x; i++){
    for(int j=0; j<n_y; j++){
      double mu_x = x.col(i).sum() / n_rows;
      double mu_y = y.col(j).sum() / n_rows;
      double s=0;

      for(int k=0; k<n_rows; k++){
        s += (x(k, i) - mu_x) * (y(k, j) - mu_y);
      }

      out_mat(i, j) = s / (n_rows - 1);
    }
  }
  return(out_mat);
}

// // [[Rcpp::export]]
// NumericMatrix normalise_cpp(NumericMatrix x){
//   int n_cols = x.cols();
//   int n_rows = x.rows();
//
//   NumericMatrix out_mat(n_rows, n_cols);
//   double mu_x = sum(x) / (n_cols*n_rows);
//   double s = 0;
//   for(int i=0; i<n_rows; i++){
//     for(int j=0; j<n_cols; j++){
//       s += (x(i, j) - mu_x) * (x(i, j) - mu_x);
//     }
//   }
//   s /= (n_cols*n_rows) - 1;
//   s = sqrt(s);
//
//   for(int i=0; i<n_rows; i++){
//     for(int j=0; j<n_cols; j++){
//       out_mat(i, j) = (x(i, j) - mu_x) / s;
//     }
//   }
//
//   return(out_mat);
// }

// [[Rcpp::export]]
Eigen::MatrixXd normalise_cpp2(Eigen::MatrixXd x){
  int n_cols = x.cols();
  int n_rows = x.rows();

  Eigen::MatrixXd out_mat(n_rows, n_cols);

  for(int i=0; i<n_cols; i++){
    double mu_x = x.col(i).sum() / n_rows;
    double s = 0;
    for(int j=0; j<n_rows; j++){
      s += (x(j, i) - mu_x) * (x(j, i) - mu_x);
    }
    s /= n_rows - 1;
    s = sqrt(s);

    for(int j=0; j<n_rows; j++){
      out_mat(j, i) = (x(j, i) - mu_x) / s;
    }
  }



  return(out_mat);
}


// [[Rcpp::export]]
Eigen::MatrixXd make_rff_cpp(Eigen::MatrixXd x, Eigen::MatrixXd w,
                             double sigma, Eigen::VectorXd b){
  w =  w / sigma;
  Eigen::MatrixXd wx = w * x.transpose();
  for(int i=0; i<wx.rows(); i++){
    for(int j=0; j<wx.cols(); j++){
      wx(i, j) += b(i);
    }
  }
  Eigen::MatrixXd outmat(wx.cols(), wx.rows());
  for(int i=0; i<wx.rows(); i++){
    for(int j=0; j<wx.cols(); j++){
      outmat(j, i) = sqrt(2) * cos(wx(i, j));
    }
  }
  return(outmat);
}

double makeSta(Eigen::MatrixXd x, Eigen::MatrixXd y,
               double median_x, double median_y,
               Eigen::MatrixXd w, Eigen::VectorXd b){
  x = normalise_cpp2(x);
  y = normalise_cpp2(y);

  int r=x.rows();

  Eigen::MatrixXd f_x = make_rff_cpp(x, w, median_x, b);
  Eigen::MatrixXd f_y = make_rff_cpp(y, w, median_y, b);

  f_x = normalise_cpp2(f_x);
  f_y = normalise_cpp2(f_y);

  Eigen::MatrixXd Cxy = cov_cpp2(f_x, f_y);
  // double Sta = r * Cxy.pow(2).sum();
  double Sta=0;

  for(int i=0; i < Cxy.rows(); i++){
    for(int j=0; j < Cxy.cols(); j++){
      Sta += Cxy(i, j) * Cxy(i, j);
    }
  }

  return(r*Sta);
}


// [[Rcpp::export]]
double RIT_cpp(Eigen::MatrixXd x, Eigen::MatrixXd y,
               double median_x, double median_y,
               Eigen::MatrixXd w, Eigen::VectorXd b,
               bool return_ts=0,
               double n_bs=100){
  double Sta = makeSta(x, y, median_x, median_y, w, b);

  if(return_ts) return(Sta);

  Eigen::MatrixXd x_copy;
  Eigen::PermutationMatrix<Eigen::Dynamic,Eigen::Dynamic> perm(x.rows());
  perm.setIdentity();
  double Sta_i;
  int n_above=0;

  for(int i=0; i<n_bs; i++){
    std::random_shuffle(perm.indices().data(), perm.indices().data()+perm.indices().size());
    x_copy = perm * x;
    Sta_i = makeSta(x_copy, y, median_x, median_y, w, b);

    if(Sta > Sta_i) n_above++;
  }

  return(1 - (n_above / n_bs));
}

double sum_squares(Eigen::MatrixXd X){
  double s=0;
  for(int i=0; i<X.rows(); i++){
    for(int j=0; j<X.cols(); j++){
      s += X(i, j) * X(i, j);
    }
  }
  return(s);
}

// [[Rcpp::export]]
double RCIT_cpp(Eigen::MatrixXd x, Eigen::MatrixXd y, Eigen::MatrixXd z,
                double median_x, double median_y, double median_z,
                Eigen::MatrixXd w, Eigen::MatrixXd wy, Eigen::MatrixXd wz,
                Eigen::VectorXd b, Eigen::VectorXd by, Eigen::VectorXd bz,
                bool return_ts=0,
                double n_bs=100){

  Eigen::MatrixXd f_x = make_rff_cpp(x, w, median_x, b);
  Eigen::MatrixXd f_y = make_rff_cpp(y, wy, median_y, by);
  Eigen::MatrixXd f_z = make_rff_cpp(z, wz, median_z, bz);

  f_x = normalise_cpp2(f_x);
  f_y = normalise_cpp2(f_y);
  f_z = normalise_cpp2(f_z);

  Eigen::MatrixXd Cxy = cov_cpp2(f_x, f_y);
  Eigen::MatrixXd Czz = cov_cpp2(f_z, f_z);
  Eigen::MatrixXd Czx = cov_cpp2(f_z, f_x);
  Eigen::MatrixXd Czy = cov_cpp2(f_z, f_y);

  for(int i=0; i<Czz.cols(); i++){
    Czz(i, i) += 0.0000000001;
  }

  Eigen::MatrixXd e_x_z = f_z * Czz.llt().solve(Czx);
  Eigen::MatrixXd e_y_z = f_z * Czz.llt().solve(Czy);

  Eigen::MatrixXd res_x = f_x - e_x_z;
  Eigen::MatrixXd res_y = f_y - e_y_z;

  Eigen::MatrixXd Cxy_z = cov_cpp2(res_x, res_y);
  int r=x.rows();
  double Sta = r*sum_squares(Cxy_z);

  if(return_ts) return(Sta);

  Eigen::MatrixXd res_x_copy;
  Eigen::PermutationMatrix<Eigen::Dynamic,Eigen::Dynamic> perm(x.rows());
  perm.setIdentity();
  double Sta_i;
  int n_above=0;

  for(int i=0; i<n_bs; i++){
    std::random_shuffle(perm.indices().data(), perm.indices().data()+perm.indices().size());
    res_x_copy = perm * res_x;
    Cxy_z = cov_cpp2(res_x_copy, res_y);

    Sta_i = r*sum_squares(Cxy_z);
    if(Sta > Sta_i) n_above++;
  }

  return(1 - (n_above / n_bs));
}



// [[Rcpp::export]]
double RCIT_disag_cpp(Eigen::MatrixXd x, Eigen::MatrixXd y, Eigen::MatrixXd z,
                      double median_x, double median_y, double median_z,
                      Eigen::MatrixXd w, Eigen::MatrixXd wy, Eigen::MatrixXd wz,
                      Eigen::VectorXd b, Eigen::VectorXd by, Eigen::VectorXd bz,
                      IntegerVector polygon_start_index,
                      IntegerVector polygon_end_index,
                      bool population_supplied,
                      NumericVector population,
                      bool return_ts=0,
                      double n_bs=100){

  int n_x = x.size();

  Eigen::MatrixXd f_x = make_rff_cpp(x, w, median_x, b);
  Eigen::MatrixXd f_y = make_rff_cpp(y, wy, median_y, by);
  Eigen::MatrixXd f_z = make_rff_cpp(z, wz, median_z, bz);

  int n_rff_z = f_z.cols();
  int n_rff_y = f_y.cols();

  Eigen::MatrixXd f_z_polygon = Eigen::MatrixXd::Constant(n_x, n_rff_z, 0);

  double polygon_pop;
  for(int i=0; i<n_x; i++){
    int i_count=0;
    polygon_pop=0;
    for(int j=polygon_start_index(i); j<polygon_end_index(i); j++){
      if(population_supplied){
        f_z_polygon.row(i) += f_z.row(j) * population(j);
        polygon_pop += population(j);
      }else{
        f_z_polygon.row(i) += f_z.row(j);
      }
      i_count++;
    }
    if(population_supplied){
      f_z_polygon.row(i) = f_z_polygon.row(i) / polygon_pop;
    }else{
      f_z_polygon.row(i) = f_z_polygon.row(i) / (double) i_count;
    }
  }

  f_x = normalise_cpp2(f_x);
  f_y = normalise_cpp2(f_y);
  f_z = normalise_cpp2(f_z);
  f_z_polygon = normalise_cpp2(f_z_polygon);


  Eigen::MatrixXd Czz = cov_cpp2(f_z, f_z);
  Eigen::MatrixXd Czz_polygon = cov_cpp2(f_z_polygon, f_z_polygon);
  Eigen::MatrixXd Czx = cov_cpp2(f_z_polygon, f_x);
  Eigen::MatrixXd Czy = cov_cpp2(f_z, f_y);

  for(int i=0; i<Czz.cols(); i++){
    Czz(i, i) += 0.000001;
    Czz_polygon(i, i) += 0.000001;

  }

  Eigen::MatrixXd e_x_z = f_z_polygon * Czz_polygon.llt().solve(Czx);
  Eigen::MatrixXd e_y_z = f_z * Czz.llt().solve(Czy);

  Eigen::MatrixXd res_x = f_x - e_x_z;
  Eigen::MatrixXd res_y = f_y - e_y_z;

  //aggregate y residuals
  Eigen::MatrixXd res_y_polygon = Eigen::MatrixXd::Constant(n_x, n_rff_y, 0);


  for(int i=0; i<n_x; i++){
    int i_count=0;
    polygon_pop=0;
    for(int j=polygon_start_index(i); j<polygon_end_index(i); j++){
      if(population_supplied){
        res_y_polygon.row(i) += res_y.row(j) * population(j);
        polygon_pop += population(j);
      }else{
        res_y_polygon.row(i) += res_y.row(j);
      }
      i_count++;
    }
    if(population_supplied){
      res_y_polygon.row(i) = res_y_polygon.row(i) / polygon_pop;
    }else{
      res_y_polygon.row(i) = res_y_polygon.row(i) / (double) i_count;
    }
  }


  Eigen::MatrixXd Cxy_z = cov_cpp2(res_x, res_y_polygon);
  int r=x.rows();
  double Sta = sum_squares(Cxy_z);

  if(return_ts) return(Sta);

  Eigen::MatrixXd res_x_copy;
  Eigen::PermutationMatrix<Eigen::Dynamic,Eigen::Dynamic> perm(x.rows());
  perm.setIdentity();
  double Sta_i;
  int n_above=0;

  for(int i=0; i<n_bs; i++){
    std::random_shuffle(perm.indices().data(), perm.indices().data()+perm.indices().size());
    res_x_copy = perm * res_x;
    Cxy_z = cov_cpp2(res_x_copy, res_y_polygon);

    Sta_i = sum_squares(Cxy_z);
    if(Sta > Sta_i) n_above++;
  }

  return(1 - (n_above / n_bs));
}

// [[Rcpp::export]]
double RIT_disag_cpp(Eigen::MatrixXd x, Eigen::MatrixXd y,
                     double median_x, double median_y,
                     Eigen::MatrixXd w, Eigen::VectorXd b,
                     IntegerVector polygon_start_index,
                     IntegerVector polygon_end_index,
                     bool population_supplied,
                     NumericVector population,
                     bool return_ts=0,
                     double n_bs=100){

  int n_x = x.size();

  Eigen::MatrixXd f_x = make_rff_cpp(x, w, median_x, b);
  Eigen::MatrixXd f_y = make_rff_cpp(y, w, median_y, b);

  int n_rff_y = f_y.cols();

  f_x = normalise_cpp2(f_x);
  f_y = normalise_cpp2(f_y);

  //aggregate y residuals
  // Eigen::MatrixXd res_y_polygon(n_x, n_rff_y);
  // res_y_polygon.Zero();

  Eigen::MatrixXd f_y_polygon = Eigen::MatrixXd::Constant(n_x, n_rff_y, 0);


  double polygon_pop;
  for(int i=0; i<n_x; i++){
    int i_count=0;
    polygon_pop=0;
    for(int j=polygon_start_index(i); j<polygon_end_index(i); j++){
      if(population_supplied){
        f_y_polygon.row(i) += f_y.row(j) * population(j);
        polygon_pop += population(j);
      }else{
        f_y_polygon.row(i) += f_y.row(j);
      }
      i_count++;
    }
    if(population_supplied){
      f_y_polygon.row(i) = f_y_polygon.row(i) / polygon_pop;
    }else{
      f_y_polygon.row(i) = f_y_polygon.row(i) / (double) i_count;
    }
  }


  Eigen::MatrixXd Cxy_z = cov_cpp2(f_x, f_y_polygon);
  int r=x.rows();
  double Sta = sum_squares(Cxy_z);

  if(return_ts) return(r*Sta);

  Eigen::MatrixXd f_x_copy;
  Eigen::PermutationMatrix<Eigen::Dynamic,Eigen::Dynamic> perm(x.rows());
  perm.setIdentity();
  double Sta_i;
  int n_above=0;

  for(int i=0; i<n_bs; i++){
    std::random_shuffle(perm.indices().data(), perm.indices().data()+perm.indices().size());
    f_x_copy = perm * f_x;
    Cxy_z = cov_cpp2(f_x_copy, f_y_polygon);

    Sta_i = sum_squares(Cxy_z);
    if(Sta > Sta_i) n_above++;
  }

  return(1 - (n_above / n_bs));
}



// [[Rcpp::export]]
double RIT_disag_cpp_v2(Eigen::MatrixXd x, Eigen::MatrixXd y,
                        double median_x, double median_y,
                        Eigen::MatrixXd w_x, Eigen::VectorXd b_x,
                        Eigen::MatrixXd w_y, Eigen::VectorXd b_y,
                        IntegerVector polygon_start_index,
                        IntegerVector polygon_end_index,
                        NumericVector population,
                        int n_obs,
                        bool return_ts=0,
                        double n_bs=100){

  // int n_x = x.size();

  Eigen::MatrixXd f_x = make_rff_cpp(x, w_x, median_x, b_x);
  Eigen::MatrixXd f_y = make_rff_cpp(y, w_y, median_y, b_y);



  int n_rff_y = f_y.cols();

  f_x = normalise_cpp2(f_x);
  f_y = normalise_cpp2(f_y);



  //aggregate y residuals
  // Eigen::MatrixXd res_y_polygon(n_x, n_rff_y);
  // res_y_polygon.Zero();

  Eigen::MatrixXd f_y_polygon = Eigen::MatrixXd::Constant(n_obs, n_rff_y, 0);
  Eigen::MatrixXd f_x_polygon = Eigen::MatrixXd::Constant(n_obs, n_rff_y, 0);


  double polygon_pop;
  for(int i=0; i<n_obs; i++){
    int i_count=0;
    polygon_pop=0;
    for(int j=polygon_start_index(i); j<polygon_end_index(i); j++){
      f_x_polygon.row(i) += f_x.row(j) * population(j);
      polygon_pop += population(j);
      f_y_polygon.row(i) += f_y.row(j) * population(j);
    }

    f_y_polygon.row(i) = f_y_polygon.row(i) / polygon_pop;
    f_x_polygon.row(i) = f_x_polygon.row(i) / polygon_pop;
  }

  // std::cout << f_x_polygon.sum() << "\n";

  Eigen::MatrixXd Cxy_z = cov_cpp2(f_x_polygon, f_y_polygon);
  int r=x.rows();
  double Sta = sum_squares(Cxy_z);

  if(return_ts) return(r*Sta);
  Eigen::MatrixXd f_x_copy;
  Eigen::PermutationMatrix<Eigen::Dynamic,Eigen::Dynamic> perm(f_x_polygon.rows());
  perm.setIdentity();
  double Sta_i;
  int n_above=0;

  for(int i=0; i<n_bs; i++){
    std::random_shuffle(perm.indices().data(), perm.indices().data()+perm.indices().size());
    f_x_copy = perm * f_x_polygon;
    // std::cout<<f_x_copy(0, 0)<<"\n";
    Cxy_z = cov_cpp2(f_x_copy, f_y_polygon);

    Sta_i = sum_squares(Cxy_z);
    // std::cout<<Sta_i<<"\n";
    if(Sta > Sta_i) n_above++;
  }

  return(1 - (n_above / n_bs));
}



// [[Rcpp::export]]
double RCIT_disag_cpp_v2(Eigen::MatrixXd x, Eigen::MatrixXd y, Eigen::MatrixXd z,
                         double median_x, double median_y, double median_z,
                         Eigen::MatrixXd w, Eigen::MatrixXd wy, Eigen::MatrixXd wz,
                         Eigen::VectorXd b, Eigen::VectorXd by, Eigen::VectorXd bz,
                         IntegerVector polygon_start_index,
                         IntegerVector polygon_end_index,
                         NumericVector population,
                         int n_obs,
                         bool return_ts=0,
                         double n_bs=100){


  Eigen::MatrixXd f_x = make_rff_cpp(x, w, median_x, b);
  Eigen::MatrixXd f_y = make_rff_cpp(y, wy, median_y, by);
  Eigen::MatrixXd f_z = make_rff_cpp(z, wz, median_z, bz);

  int n_rff_z = f_z.cols();
  int n_rff_y = f_y.cols();

  Eigen::MatrixXd f_z_polygon = Eigen::MatrixXd::Constant(n_obs, n_rff_z, 0);
  Eigen::MatrixXd f_x_polygon = Eigen::MatrixXd::Constant(n_obs, n_rff_y, 0);
  Eigen::MatrixXd f_y_polygon = Eigen::MatrixXd::Constant(n_obs, n_rff_y, 0);

  double polygon_pop;
  for(int i=0; i<n_obs; i++){
    int i_count=0;
    polygon_pop=0;
    for(int j=polygon_start_index(i); j<polygon_end_index(i); j++){
      f_z_polygon.row(i) += f_z.row(j) * population(j);
      f_x_polygon.row(i) += f_x.row(j) * population(j);
      f_y_polygon.row(i) += f_y.row(j) * population(j);
      polygon_pop += population(j);
    }
    f_z_polygon.row(i) = f_z_polygon.row(i) / polygon_pop;
    f_x_polygon.row(i) = f_x_polygon.row(i) / polygon_pop;
    f_y_polygon.row(i) = f_y_polygon.row(i) / polygon_pop;
  }

  f_x_polygon = normalise_cpp2(f_x_polygon);
  f_y_polygon = normalise_cpp2(f_y_polygon);
  f_z_polygon = normalise_cpp2(f_z_polygon);


  Eigen::MatrixXd Czz = cov_cpp2(f_z_polygon, f_z_polygon);
  Eigen::MatrixXd Czx = cov_cpp2(f_z_polygon, f_x_polygon);
  Eigen::MatrixXd Czy = cov_cpp2(f_z_polygon, f_y_polygon);

  for(int i=0; i<Czz.cols(); i++){
    Czz(i, i) += 0.000001;
  }

  Eigen::MatrixXd e_x_z = f_z_polygon * Czz.llt().solve(Czx);
  Eigen::MatrixXd e_y_z = f_z_polygon * Czz.llt().solve(Czy);

  Eigen::MatrixXd res_x = f_x - e_x_z;
  Eigen::MatrixXd res_y = f_y - e_y_z;

  Eigen::MatrixXd Cxy_z = cov_cpp2(res_x, res_y);
  int r=x.rows();
  double Sta = sum_squares(Cxy_z);

  if(return_ts) return(r*Sta);

  Eigen::MatrixXd res_x_copy;
  Eigen::PermutationMatrix<Eigen::Dynamic,Eigen::Dynamic> perm(n_obs);
  perm.setIdentity();
  double Sta_i;
  int n_above=0;

  for(int i=0; i<n_bs; i++){
    std::random_shuffle(perm.indices().data(), perm.indices().data()+perm.indices().size());
    res_x_copy = perm * res_x;
    Cxy_z = cov_cpp2(res_x_copy, res_y);

    Sta_i = sum_squares(Cxy_z);
    if(Sta > Sta_i) n_above++;
  }

  return(1 - (n_above / n_bs));
}
