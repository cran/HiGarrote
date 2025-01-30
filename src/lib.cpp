// [[Rcpp::plugins("cpp11")]]
#include <RcppArmadillo.h>

//////////////////////////////// get_level
// [[Rcpp::export]]
Rcpp::List get_level(Rcpp::DataFrame D, int p) {
  Rcpp::List uni_list(p);
  for (int i = 0; i < p; ++i) {
    Rcpp::NumericVector v = D[i];
    Rcpp::NumericVector uni_v = Rcpp::unique(v);
    std::sort(uni_v.begin(), uni_v.end());
    uni_list[i] = uni_v;
  }
  return uni_list;
}


//////////////////////////////// contr_scale
// [[Rcpp::export]]
Rcpp::NumericMatrix contr_scale(Rcpp::NumericMatrix x, int level_num) {
  arma::mat x_arma(x.begin(), x.nrow(), x.ncol(), false);
  arma::mat x2 = x_arma.t() * x_arma;
  arma::vec x3 = x2.diag()/level_num;
  arma::vec x4 = 1.0 / arma::sqrt(x3);
  arma::mat x5 = x_arma * (arma::diagmat(x4));
  Rcpp::NumericMatrix x6 = Rcpp::wrap(x5);
  return x6;
}


//////////////////////////////// U_j_cpp
// [[Rcpp::export]]
Rcpp::List U_j_cpp(Rcpp::List uni_level, int p, Rcpp::IntegerVector mi,
                   SEXP quali_id, SEXP quanti_eq_id, SEXP quanti_ineq_id,
                   SEXP quali_contr) {
  
  Rcpp::List U_j_list(p);
  
  Rcpp::Environment statsEnv = Rcpp::Environment::namespace_env("stats");
  Rcpp::Function poly = statsEnv["poly"];
  Rcpp::Function contrasts = statsEnv["contrasts"];
  Rcpp::Function contr_helmert = statsEnv["contr.helmert"];
  Rcpp::Function contr_poly = statsEnv["contr.poly"];
  
  Rcpp::IntegerVector quali_id_cpp;
  Rcpp::IntegerVector quanti_eq_id_cpp;
  Rcpp::IntegerVector quanti_ineq_id_cpp;
  Rcpp::List quali_contr_cpp;
  
  if(!Rf_isNull(quali_id)) {
    quali_id_cpp = Rcpp::clone(Rcpp::as<Rcpp::IntegerVector>(quali_id));
    quali_id_cpp = quali_id_cpp - 1;
  } 
  if(!Rf_isNull(quanti_eq_id)) {
    quanti_eq_id_cpp = Rcpp::clone(Rcpp::as<Rcpp::IntegerVector>(quanti_eq_id));
    quanti_eq_id_cpp = quanti_eq_id_cpp - 1;
  }  
  if(!Rf_isNull(quanti_ineq_id)) {
    quanti_ineq_id_cpp = Rcpp::clone(Rcpp::as<Rcpp::IntegerVector>(quanti_ineq_id));
    quanti_ineq_id_cpp = quanti_ineq_id_cpp - 1;
  }
  if(!Rf_isNull(quali_contr)) {
    quali_contr_cpp = Rcpp::clone(Rcpp::as<Rcpp::List>(quali_contr));
  }
  
  for (int j = 0; j < p; ++j) {
    Rcpp::NumericVector x = uni_level[j];
    Rcpp::NumericMatrix x2;
    
    if (std::find(quali_id_cpp.begin(), quali_id_cpp.end(), j) != quali_id_cpp.end()) {
      Rcpp::NumericMatrix x1;
      if(!Rf_isNull(quali_contr_cpp[j])) {
        x1 = Rcpp::as<Rcpp::NumericMatrix>(quali_contr_cpp[j]);
      } else {
        x1 = contr_helmert(x);
      }
      x2 = contr_scale(x1, mi[j]);
      arma::mat x_arma(x2.begin(), x2.nrow(), x2.ncol(), false);
      int n_rows = x_arma.n_rows;
      arma::vec ones_column = arma::ones<arma::vec>(n_rows);
      arma::mat new_matrix = arma::join_horiz(ones_column, x_arma);
      x2 = Rcpp::wrap(new_matrix);
      
    } else if (std::find(quanti_eq_id_cpp.begin(), quanti_eq_id_cpp.end(), j) != quanti_eq_id_cpp.end()) {
      Rcpp::NumericMatrix x1 = contr_poly(x);
      x2 = contr_scale(x1, mi[j]);
      arma::mat x_arma(x2.begin(), x2.nrow(), x2.ncol(), false);
      int n_rows = x_arma.n_rows;
      arma::vec ones_column = arma::ones<arma::vec>(n_rows);
      arma::mat new_matrix = arma::join_horiz(ones_column, x_arma);
      x2 = Rcpp::wrap(new_matrix);
      
    } else if (std::find(quanti_ineq_id_cpp.begin(), quanti_ineq_id_cpp.end(), j) != quanti_ineq_id_cpp.end()) {
      Rcpp::NumericMatrix x1 = poly(x, mi[j] - 1);
      x2 = contr_scale(x1, mi[j]);
      arma::mat x_arma(x2.begin(), x2.nrow(), x2.ncol(), false);
      int n_rows = x_arma.n_rows;
      arma::vec ones_column = arma::ones<arma::vec>(n_rows);
      arma::mat new_matrix = arma::join_horiz(ones_column, x_arma);
      x2 = Rcpp::wrap(new_matrix);
      
    } else {
      arma::mat x1(2,2);
      x1(0,0) = 1;
      x1(1,0) = 1;
      x1(0,1) = -1;
      x1(1,1) = 1;
      x2 = Rcpp::wrap(x1);
    }
    
    
    U_j_list[j] = x2;
    
    // Rcpp::CharacterMatrix m_char = Rcpp::as<Rcpp::CharacterMatrix>(m);
    // U_j_list[j] = m_char;
  }
  return U_j_list;
}


//////////////////////////////// arma_dist
arma::mat arma_dist(const arma::vec& x) {
  int n = x.n_elem;
  arma::mat dist_mat(n,n,arma::fill::zeros);
  
  for(int i = 0; i < n; ++i) {
    for(int j = i+1; j < n; ++j) {
      double d = std::abs(x(i) - x(j));
      dist_mat(i,j) = d;
      dist_mat(j,i) = d;
    }
  }
  return dist_mat;
}


//////////////////////////////// h_dist_cpp
// [[Rcpp::export]]
Rcpp::List h_dist_cpp(Rcpp::NumericVector x, Rcpp::NumericMatrix m,
                      bool two_level, bool quali) {
  arma::vec x_arma = Rcpp::as<arma::vec>(x);

  if(two_level) {
    Rcpp::List h_list(1);
    arma::mat h = arma_dist(x);
    h.transform([](double val) {return val != 0 ? 1.0 : 0.0;});
    h_list[0] = h;
    return Rcpp::List::create(h_list);

  } else if(quali) {
    arma::mat m_arma(m.begin(), m.nrow(), m.ncol(), false);
    arma::mat m_arma1 = m_arma.cols(1, m_arma.n_cols-1);
    
    int k = m_arma1.n_cols;
    Rcpp::List h_list(k);
    for (int i = 0; i < k; ++i) {
      arma::vec col = m_arma1.col(i);
      arma::mat h_col = arma_dist(col);
      h_col.transform( [](double val) { return val != 0 ? 1.0 : 0.0; } );
      h_list[i] = h_col;
    }
    return Rcpp::List::create(h_list);
    
  } else {
    Rcpp::List h_list(1);
    arma::mat h = arma_dist(x_arma);
    h_list[0] = h;
    return Rcpp::List::create(h_list);
  }
}

//////////////////////////////// h_j_cpp
// [[Rcpp::export]]
Rcpp::List h_j_cpp(int p, Rcpp::List uni_level, Rcpp::List m_list,
                   SEXP two_level_id, SEXP quali_id) {
  Rcpp::List h_j_list(p);
  
  Rcpp::NumericVector two_level_id_cpp;
  Rcpp::NumericVector quali_id_cpp;
  if(!Rf_isNull(two_level_id)) {
    two_level_id_cpp = Rcpp::clone(Rcpp::as<Rcpp::NumericVector>(two_level_id));
    two_level_id_cpp = two_level_id_cpp - 1;
  }
  if(!Rf_isNull(quali_id)) {
    quali_id_cpp = Rcpp::clone(Rcpp::as<Rcpp::NumericVector>(quali_id));
    quali_id_cpp = quali_id_cpp - 1;
  }
  
  for(int j = 0; j < p; ++j) {
    if(std::find(two_level_id_cpp.begin(), two_level_id_cpp.end(), j) != two_level_id_cpp.end()) {
      h_j_list[j] = h_dist_cpp(uni_level[j], m_list[j], true, false);
    } else if(std::find(quali_id_cpp.begin(), quali_id_cpp.end(), j) != quali_id_cpp.end()) {
      h_j_list[j] = h_dist_cpp(uni_level[j], m_list[j], false, true);
    } else {
      h_j_list[j] = h_dist_cpp(uni_level[j], m_list[j], false, false);
    }
  }
  return h_j_list;
}

//////////////////////////////// Psi_mat_cpp
// [[Rcpp::export]]
arma::mat Psi_mat_cpp(const std::vector<arma::mat>& h_list_mat, const Rcpp::NumericVector& rho) {
  arma::mat first_mat = h_list_mat[0];
  int n = first_mat.n_rows;
  arma::mat Psi = arma::ones(n, n);
  
  for(int i = 0; i < h_list_mat.size(); ++i) {
    arma::mat h_list_mat_i = h_list_mat[i];
    double rho_i = rho[i];
    Psi %= arma::exp(std::log(rho_i) * arma::square(h_list_mat_i)); 
  }
  return Psi;
}


//////////////////////////////// NLLH
class NLLH {
public:
  NLLH(const Rcpp::List& h_list_mat_, int n_, int replicate_, const Rcpp::NumericVector& y_,
       double nugget_ = 1e-6, double epsilon_ = 1e-3, bool interpolate_ = true)
    : n(n_), replicate(replicate_), y(y_), nugget(nugget_), epsilon(epsilon_), interpolate(interpolate_){
    
    for(int i = 0; i < h_list_mat_.size(); ++i) {
      arma::mat mat = Rcpp::as<arma::mat>(h_list_mat_[i]);
      h_list_mat.push_back(mat);
    }
    
    y = Rcpp::as<arma::vec>(y_);
  }
  
  // Function to get a list excluding the i-th matrix using iterators
  std::vector<arma::mat> exclude_i_mat(int i) {
    std::vector<arma::mat> result;
    result.reserve(h_list_mat.size() - 1);
    result.insert(result.end(), h_list_mat.begin(), h_list_mat.begin() + i);
    result.insert(result.end(), h_list_mat.begin() + i + 1, h_list_mat.end());
    return result;
  }
  
  Rcpp::NumericVector exclude_i_ele(const Rcpp::NumericVector& rho, int i) {
    Rcpp::NumericVector result(rho.size() - 1);
    for (int j = 0, k = 0; j < rho.size(); ++j) {
      if (j != i) {
        result[k++] = rho[j];  // Only copy elements that are not the i-th element
      }
    }
    return result;
  }
  
  //////////////////////////////// nllh_cpp
  Rcpp::List nllh_cpp(const Rcpp::NumericVector& rho_lambda) {
    Rcpp::NumericVector rho = rho_lambda[Rcpp::Range(0, rho_lambda.size()-2)];
    int rho_size = rho.size();
    double lambda = rho_lambda[rho_lambda.size() - 1];
    
    lambda = std::min(std::max(lambda, 1e-12), 0.999);
    // Ensure valid rho and lambda values using Rcpp::all() and Rcpp::is_finite()
    if (!Rcpp::all(Rcpp::is_finite(rho)).is_true() || !std::isfinite(lambda)) {
      // Return Inf for invalid objective and gradient to signal a bad point in optimization
      return Rcpp::List::create(Rcpp::Named("objective") = R_PosInf, Rcpp::Named("gradient") = Rcpp::NumericVector(rho_size + 1, R_PosInf));
    }
     
    arma::mat Psi_mat = Psi_mat_cpp(h_list_mat, rho);
    arma::mat R_beta = Psi_mat + (lambda/(1-lambda))*(arma::eye(n,n)/replicate);
    
    // Eigen decomposition
    arma::vec eigval;
    arma::mat eigvec;
    arma::eig_sym(eigval, eigvec, R_beta);

    // Cholesky decomposition
    arma::vec a = eigvec.t() * arma::ones(n);
    arma::vec b = eigvec.t() * y;
    double mu = arma::sum(a%b / eigval) / arma::sum(arma::square(a) / eigval);
    double nu2 = (1.0/n) * (arma::sum(arma::square(b) / eigval) - mu*mu*arma::sum(arma::square(a) / eigval));
    nu2 = std::max(nu2, 1e-15);

    // obj
    double obj = std::log(nu2) + arma::sum(arma::log(eigval))/n;

    // gradient
    arma::mat R_beta_inv = eigvec * arma::diagmat(1/eigval) * eigvec.t();

    std::vector<arma::mat> R_beta_deriv(rho_size);
    std::vector<arma::mat> h_list_mat_no_i;
    Rcpp::NumericVector rho_no_i;
    for(int i = 0; i < rho_size; ++i) {
      arma::mat g1 = arma::square(h_list_mat[i]) % arma::exp(std::log(rho[i]) * (arma::square(h_list_mat[i])-1));
      h_list_mat_no_i = exclude_i_mat(i);
      rho_no_i = exclude_i_ele(rho, i);
      arma::mat g2 = Psi_mat_cpp(h_list_mat_no_i, rho_no_i);
      R_beta_deriv[i] = g1 % g2;
    }
    arma::mat lambda_deriv = (1.0/std::pow(1.0-lambda, 2)) * (arma::eye(n,n)/replicate);
    R_beta_deriv.push_back(lambda_deriv);

    arma::vec gradient(rho_size + 1);
    for(int i = 0; i < R_beta_deriv.size(); ++i) {
      double part1 = arma::as_scalar(y.t()*R_beta_inv*R_beta_deriv[i]*R_beta_inv*y) / arma::as_scalar(y.t()*R_beta_inv*y);
      double part2 = arma::trace(R_beta_inv*R_beta_deriv[i])/n;
      gradient(i) = -(part1 - part2);
    }

    return Rcpp::List::create(Rcpp::Named("objective") = obj, Rcpp::Named("gradient") = gradient);
  }
  
  // //////////////////////////////// small_GP_cpp
  // Rcpp::List small_GP_cpp(const Rcpp::NumericVector& rho_lambda) {
  //   Rcpp::NumericVector rho = rho_lambda[Rcpp::Range(0, rho_lambda.size()-2)];
  //   double lambda = rho_lambda[rho_lambda.size() - 1];
  //   
  //   arma::mat Psi_mat = Psi_mat_cpp(h_list_mat, rho);
  //   arma::mat R_beta0 = Psi_mat + (lambda/(1-lambda))*(arma::eye(n,n));
  //   arma::mat R_beta = Psi_mat + (lambda/(1-lambda))*(arma::eye(n,n)/replicate);
  //   
  //   // Eigen decomposition
  //   arma::vec eigval;
  //   arma::mat eigvec;
  //   arma::eig_sym(eigval, eigvec, R_beta);
  //   
  //   // Cholesky decomposition
  //   arma::vec a = eigvec.t() * arma::ones(n);
  //   arma::vec b = eigvec.t() * y;
  //   double mu = arma::sum(a%b / eigval) / arma::sum(arma::square(a) / eigval);
  //   double nu2 = (1.0/n) * (arma::sum(arma::square(b) / eigval) - mu*mu*arma::sum(arma::square(a) / eigval));
  //   nu2 = std::max(nu2, 1e-15);
  //   
  //   return Rcpp::List::create(Rcpp::Named("mu") = mu, Rcpp::Named("nu2") = nu2,
  //                             Rcpp::Named("a") = a, Rcpp::Named("b") = b, 
  //                             Rcpp::Named("eigval") = eigval, Rcpp::Named("eigvec") = eigvec,
  //                             Rcpp::Named("Psi_mat") = Psi_mat, Rcpp::Named("R_beta") = R_beta);
  // }
  
  //////////////////////////////// nllh_GP_cpp
  Rcpp::List nllh_GP_cpp(const Rcpp::NumericVector& rho) {
    int rho_size = rho.size();
    
    // Ensure valid rho values using Rcpp::all() and Rcpp::is_finite()
    if (!Rcpp::all(Rcpp::is_finite(rho)).is_true() ) {
      // Return Inf for invalid objective and gradient to signal a bad point in optimization
      return Rcpp::List::create(Rcpp::Named("objective") = R_PosInf, Rcpp::Named("gradient") = Rcpp::NumericVector(rho_size, R_PosInf));
    }
    
    arma::mat Psi_mat = Psi_mat_cpp(h_list_mat, rho);
    arma::mat R_beta = Psi_mat + nugget*(arma::eye(n,n)/replicate);
    
    // Eigen decomposition
    arma::vec eigval;
    arma::mat eigvec;
    arma::eig_sym(eigval, eigvec, R_beta);
    
    // Cholesky decomposition
    arma::vec a = eigvec.t() * arma::ones(n);
    arma::vec b = eigvec.t() * y;
    arma::mat L = arma::chol(R_beta, "lower");
    arma::vec a0 = arma::solve(arma::trimatl(L), arma::ones(n), arma::solve_opts::fast);
    arma::vec b0 = arma::solve(arma::trimatl(L), y, arma::solve_opts::fast);
    // arma::mat L = eigvec * arma::diagmat(arma::sqrt(eigval));
    // arma::vec a0 = a * (1/arma::sqrt(eigval));
    // arma::vec b0 = b * (1/arma::sqrt(eigval));
    
    double mu = arma::sum(a%b / eigval) / arma::sum(arma::square(a) / eigval);
    double nu2 = (1.0/(n-1)) * (arma::sum(arma::square(b) / eigval) - mu*mu*arma::sum(arma::square(a) / eigval));
    nu2 = std::max(nu2, 1e-15);
    
    // obj
    // double obj = (n-1)*std::log(nu2) + arma::sum(arma::log(eigval)) + std::log(arma::sum(arma::square(a0)));
    double gg1 = (n-1)*std::log(nu2);
    double gg2 = arma::sum(arma::log(eigval));
    double gg3 = std::log(arma::sum(arma::square(a0)));
    double obj = gg1 + gg2 + gg3;
    if (interpolate) {
      // penalty for non-interpolating
      arma::vec ym = b0 - mu * a0;
      arma::vec L_t_solve_ym = arma::solve(arma::trimatu(L.t()), ym, arma::solve_opts::fast);
      arma::vec pred = mu + (L*ym - L_t_solve_ym*nugget);
      double y_tss = arma::sum(arma::square(y)) - std::pow(arma::sum(y), 2)/n;
      double penalty = 2 * arma::sum(arma::square((y-pred))) / (epsilon * y_tss);
      obj = obj / n + penalty;
    }
    
    
    // gradient
    arma::mat R_beta_inv = eigvec * arma::diagmat(1/eigval) * eigvec.t();
    
    std::vector<arma::mat> R_beta_deriv(rho_size);
    std::vector<arma::mat> h_list_mat_no_i;
    Rcpp::NumericVector rho_no_i;
    for(int i = 0; i < rho_size; ++i) {
      arma::mat g1 = arma::square(h_list_mat[i]) % arma::exp(std::log(rho[i]) * (arma::square(h_list_mat[i])-1));
      h_list_mat_no_i = exclude_i_mat(i);
      rho_no_i = exclude_i_ele(rho, i);
      arma::mat g2 = Psi_mat_cpp(h_list_mat_no_i, rho_no_i);
      R_beta_deriv[i] = g1 % g2;
    }
    arma::vec gradient(rho_size);
    for(int i = 0; i < R_beta_deriv.size(); ++i) {
      double part1 = arma::as_scalar(y.t()*R_beta_inv*R_beta_deriv[i]*R_beta_inv*y) / arma::as_scalar(y.t()*R_beta_inv*y);
      double part2 = arma::trace(R_beta_inv*R_beta_deriv[i])/n;
      gradient(i) = -(part1 - part2);
    }

    return Rcpp::List::create(Rcpp::Named("objective") = obj, Rcpp::Named("gradient") = gradient,
                              Rcpp::Named("mu") = mu, Rcpp::Named("nu2") = nu2,
                              Rcpp::Named("a") = a, Rcpp::Named("b") = b,
                              Rcpp::Named("eigval") = eigval, Rcpp::Named("eigvec") = eigvec);
  } 
  
  
private:
  std::vector<arma::mat> h_list_mat;
  int n;
  int replicate;
  arma::vec y;
  double nugget;
  double epsilon;
  bool interpolate;
};

//////////////////////////////// Export nllh_cpp
// Define a global pointer to an NLLH instance
Rcpp::XPtr<NLLH>* NLLH_instance_ptr = nullptr;

// Function to initialize the NLLH instance
// [[Rcpp::export]]
void initialize_NLLH_instance(Rcpp::List h_list_mat, int n, int replicate, Rcpp::NumericVector y,
                              double nugget, double epsilon, bool interpolate) {
  NLLH* ptr = new NLLH(h_list_mat, n, replicate, y, nugget, epsilon, interpolate);
  NLLH_instance_ptr = new Rcpp::XPtr<NLLH>(ptr, true);
} 

// Function to calculate the negative log-likelihood
// [[Rcpp::export]]
Rcpp::List nllh_cpp_R(Rcpp::NumericVector rho_lambda) {
  if (NLLH_instance_ptr == nullptr) {
    Rcpp::stop("NLLH instance is not initialized. Call initialize_NLLH_instance first.");
  }
  return NLLH_instance_ptr->get()->nllh_cpp(rho_lambda);
}

// Function to get information required in GP
// [[Rcpp::export]]
Rcpp::List nllh_GP_R(Rcpp::NumericVector rho) {
  if (NLLH_instance_ptr == nullptr) {
    Rcpp::stop("NLLH instance is not initialized. Call initialize_NLLH_instance first.");
  }
  return NLLH_instance_ptr->get()->nllh_GP_cpp(rho);
}


//////////////////////////////// rho_lambda_optim
// [[Rcpp::export]]
Rcpp::List rho_lambda_optim(const Rcpp::NumericMatrix& ini_point, const Rcpp::List& h_list_mat,
                            int n, int replicate, Rcpp::NumericVector y,
                            double lambda_lb, double lambda_ub, 
                            double nugget = 1e-6, double epsilon = 1e-3, bool interpolate = true) {
  initialize_NLLH_instance(h_list_mat, n, replicate, y, nugget, epsilon, interpolate);
  int num_point = ini_point.nrow();
  int num_param = ini_point.ncol();
  Rcpp::List result(num_point);

  Rcpp::Environment nloptrEnv = Rcpp::Environment::namespace_env("nloptr");
  Rcpp::Function nloptr = nloptrEnv["nloptr"];

  // Set up
  Rcpp::List opts = Rcpp::List::create(
    Rcpp::Named("algorithm") = "NLOPT_LD_MMA",
    Rcpp::Named("xtol_rel") = 1.0e-8,
    Rcpp::Named("maxeval") = 100
  );
  arma::vec lb = arma::join_cols(arma::vec(num_param - 1).fill(1e-15), arma::vec({lambda_lb}));
  arma::vec ub = arma::join_cols(arma::vec(num_param - 1).fill(0.999), arma::vec({lambda_ub}));

  Rcpp::NumericVector ini_point_i;
  for(int i = 0; i < num_point; ++i) {
    ini_point_i = ini_point(i,Rcpp::_);
    Rcpp::List x = nloptr(Rcpp::_["x0"] = ini_point_i,
                          Rcpp::_["eval_f"] = Rcpp::InternalFunction(&nllh_cpp_R),
                          Rcpp::_["lb"] = lb,
                          Rcpp::_["ub"] = ub,
                          Rcpp::_["opts"] = opts);
    result[i] = Rcpp::List::create(x);
  }

  return result;
}


//////////////////////////////// rho_optim_GP
// [[Rcpp::export]]
Rcpp::List rho_optim_GP(const Rcpp::NumericMatrix& ini_point, const Rcpp::List& h_list_mat,
                        int n, int replicate, Rcpp::NumericVector y,
                        double nugget = 1e-6, double epsilon = 1e-3, bool interpolate = true) {
  initialize_NLLH_instance(h_list_mat, n, replicate, y, nugget, epsilon, interpolate);
  
  int num_point = ini_point.nrow();
  int num_param = ini_point.ncol();
  Rcpp::List result(num_point);

  Rcpp::Environment nloptrEnv = Rcpp::Environment::namespace_env("nloptr");
  Rcpp::Function nloptr = nloptrEnv["nloptr"];

  // Set up
  Rcpp::List opts = Rcpp::List::create(
    Rcpp::Named("algorithm") = "NLOPT_LD_MMA",
    Rcpp::Named("xtol_rel") = 1.0e-8,
    Rcpp::Named("maxeval") = 100
  );
  arma::vec lb = arma::vec(num_param).fill(1e-15);
  arma::vec ub = arma::vec(num_param).fill(0.999);
  
  Rcpp::NumericVector ini_point_i;
  for(int i = 0; i < num_point; ++i) {
    ini_point_i = ini_point(i,Rcpp::_);
    Rcpp::List x = nloptr(Rcpp::_["x0"] = ini_point_i,
                          Rcpp::_["eval_f"] = Rcpp::InternalFunction(&nllh_GP_R),
                          Rcpp::_["lb"] = lb,
                          Rcpp::_["ub"] = ub,
                          Rcpp::_["opts"] = opts);
    result[i] = Rcpp::List::create(x);
  } 
  
  return result;
} 


//////////////////////////////// Psi_j_list_cpp
std::vector<arma::mat> Psi_j_list_cpp(const Rcpp::List& h_j_list, const Rcpp::List& rho_list, int p) {
  std::vector<arma::mat> Psi_j_list(p);
  
  for(int j = 0; j < p; ++j) {
    Rcpp::List h_j_list_mat = Rcpp::as<Rcpp::List>(h_j_list[j]);
    std::vector<arma::mat> h_j_list_mat_arma;
    
    for(int k = 0; k < h_j_list_mat.size(); ++k) {
      arma::mat a = Rcpp::as<arma::mat>(h_j_list_mat[k]);
      h_j_list_mat_arma.push_back(a);
    }
    
    // Psi_j
    Rcpp::NumericVector rho_j = rho_list[j];
    arma::mat Psi_j = Psi_mat_cpp(h_j_list_mat_arma, rho_j);
    Psi_j_list[j] = Psi_j;
  } 
  return Psi_j_list;
} 


//////////////////////////////// tau0_sigma0_cpp
double tau0_sigma0_cpp(const std::vector<arma::mat> Psi_j_list, int p, Rcpp::IntegerVector mi) {
  arma::mat Psi_j;
  arma::vec Psi_j_sum(p);
  
  for(int j = 0; j < p; ++j) {
    Psi_j = Psi_j_list[j];
    Psi_j_sum(j) = arma::accu(Psi_j);
  }
  
  double part1 = arma::prod(Psi_j_sum);
  arma::vec mi_arma = Rcpp::as<arma::vec>(mi);
  double part2 = arma::prod(mi_arma);
  double result = part1/std::pow(part2, 2);
  
  return result;
} 

// // is_PD
// // [[Rcpp::export]]
// bool is_PD(const arma::mat& M) {
//   // Attempt Cholesky decomposition
//   arma::mat L;
//   return arma::chol(L, M);  // Returns false if not positive definite
// }

// // make_PD
// arma::mat make_PD(const arma::mat& M) {
//   arma::vec eigval;
//   arma::mat eigvec;
//   arma::eig_sym(eigval, eigvec, M);
// 
//   eigval.elem(arma::find(eigval < 1e-6)).fill(1e-6);
//   arma::mat M_pd = eigvec * arma::diagmat(eigval) * eigvec.t();
//   return M_pd;
// }


//////////////////////////////// BETA
class BETA {
public:
  BETA(const Rcpp::List& h_j_list_, int p_, const Rcpp::List& rho_list_, Rcpp::IntegerVector mi_)
    : h_j_list(h_j_list_), p(p_), rho_list(rho_list_), mi(mi_) {
    
    Psi_j_list = Psi_j_list_cpp(h_j_list, rho_list, p); // Psi_j_list
    tau0_sigma0 = tau0_sigma0_cpp(Psi_j_list, p, mi); // tau0_sigma0
  }
  
  //////////////////////////////// r_j_cpp
  Rcpp::List r_j_cpp(const Rcpp::List& U_j_list, const Rcpp::IntegerVector& me_num) {
    arma::mat U_j;
    arma::mat Psi_j; 
    arma::mat a;
    arma::vec r_j_all; 
    int me_num_j;
    Rcpp::List result(p);
    
    for(int j = 0; j < p; ++j) {
      U_j = Rcpp::as<arma::mat>(U_j_list[j]);
      Psi_j = Psi_j_list[j];
      a = U_j.t() * Psi_j * U_j;
      r_j_all = a.diag();
      
      me_num_j = me_num[j];
      arma::vec b(me_num_j);
      for(int k = 0; k < me_num_j; k++) {
        b(k) = r_j_all(k+1)/r_j_all(0);
      }
      result[j] = b;
    } 
    return result;
  } 
  
  //////////////////////////////// beta_ele_cpp
  Rcpp::List beta_ele_cpp(const Rcpp::NumericMatrix& U, const Rcpp::NumericVector& R,
                          double lambda, int replicate, int n,
                          const Rcpp::NumericVector& y) {
    arma::mat U_arma = Rcpp::as<arma::mat>(U);
    arma::vec R_arma = Rcpp::as<arma::vec>(R);
    arma::vec y_arma = Rcpp::as<arma::vec>(y);
    
    arma::mat U_R_U0 = U_arma.each_row() % R_arma.t();
    arma::mat U_R_U = U_R_U0 * U_arma.t();
    arma::mat diag_term = (lambda/(1-lambda)) * (arma::eye(n,n)/replicate);
    arma::mat R_beta = tau0_sigma0*U_R_U + diag_term;
    
    // Cholesky decomposition
    arma::mat L = arma::chol(R_beta, "lower");
    arma::vec a = arma::solve(arma::trimatl(L), arma::ones<arma::vec>(n));
    arma::vec b = arma::solve(arma::trimatl(L), y_arma);
    double mu = arma::sum(a%b) / arma::sum(arma::square(a));
    double nu2 = (arma::sum(arma::square(b)) - mu*arma::sum(a%b))/n;
    nu2 = std::max(nu2, 1e-15);
    arma::vec R_beta_inv_y_mu = arma::solve(arma::trimatu(arma::trans(L)), (b-mu*a));
    arma::vec beta_arma = tau0_sigma0 * (U_R_U0.t() * R_beta_inv_y_mu);
    
    // beta_nng
    arma::mat Z = U_arma.each_row() % beta_arma.t();
    arma::mat D_nng = Z.t() * Z;
    
    arma::vec d = Z.t() * y_arma;
    Rcpp::NumericVector d_R = Rcpp::wrap(d);
    
    return Rcpp::List::create(Rcpp::Named("Dmat") = D_nng,
                              Rcpp::Named("dvec") = d,
                              Rcpp::Named("Z") = Z,
                              Rcpp::Named("U") = U,
                              Rcpp::Named("U_R_U0") = U_R_U0,
                              Rcpp::Named("R_beta") = R_beta,
                              Rcpp::Named("beta") = beta_arma);
  } 
  
  
  //////////////////////////////// beta_nng_cpp
  Rcpp::NumericVector beta_nng_cpp(const Rcpp::List& beta_ele, int replicate, int n,
                                   const Rcpp::NumericVector& y, const Rcpp::NumericMatrix& Amat,
                                   double s2) {
    Rcpp::NumericMatrix D_nng_R = beta_ele[0];
    Rcpp::NumericVector d_R = beta_ele[1];
    arma::mat Z = Rcpp::as<arma::mat>(beta_ele[2]);
    arma::mat U_arma = Rcpp::as<arma::mat>(beta_ele[3]);
    arma::mat U_R_U0 = Rcpp::as<arma::mat>(beta_ele[4]);
    arma::mat R_beta = Rcpp::as<arma::mat>(beta_ele[5]);
    arma::mat beta_arma = Rcpp::as<arma::mat>(beta_ele[6]);
    arma::vec y_arma = Rcpp::as<arma::vec>(y);
    
    // GCV
    Rcpp::Environment quadprogEnv = Rcpp::Environment::namespace_env("quadprog");
    Rcpp::Function solve_QP = quadprogEnv["solve.QP"];

    arma::vec best_coef_nng;
    arma::vec beta_nng_arma;

    if(replicate == 1) {
      double start = 0.1;
      double end = (n-1)*0.3;
      arma::vec M = arma::linspace(start, end, 100);
      arma::vec cv(100);
      int Amat_col = Amat.ncol();
      Rcpp::NumericVector b0 = Rcpp::rep(0.0, Amat_col);
      Rcpp::List coef_nng_list(100);

      for(int i = 0; i < 100; ++i) {
        b0[0] = -M(i);
        Rcpp::List result = solve_QP(Rcpp::Named("Dmat") = D_nng_R,
                                     Rcpp::Named("dvec") = d_R,
                                     Rcpp::Named("Amat") = Amat,
                                     Rcpp::Named("bvec") = b0);
        Rcpp::NumericVector coef_nng_R = result["solution"];
        arma::vec coef_nng = Rcpp::as<arma::vec>(coef_nng_R);
        arma::vec e = y_arma - Z*coef_nng;
        arma::mat df_1 = U_arma.each_row() % coef_nng.t();
        arma::mat df_2 = df_1 * U_R_U0.t() * arma::solve(R_beta, arma::eye(n,n));
        double df = tau0_sigma0 * arma::trace(df_2);
        double xx = std::pow(1-df/n, 2);
        cv(i) = arma::sum(arma::square(e)) / (n*xx);
        coef_nng_list[i] = coef_nng_R;
      }
      int best_idx = arma::index_min(cv);
      best_coef_nng = Rcpp::as<arma::vec>(coef_nng_list[best_idx]);


    } else {
      double nng_lambda = s2;
      arma::vec d = Rcpp::as<arma::vec>(d_R);
      d -= nng_lambda;
      Rcpp::NumericVector d_R2 = Rcpp::wrap(d);
      arma::mat Amat_arma = Rcpp::as<arma::mat>(Amat);
      Amat_arma.shed_col(0);
      Rcpp::NumericMatrix Amat_R = Rcpp::wrap(Amat_arma);
      Rcpp::NumericVector b0 = Rcpp::rep(0.0, Amat_R.ncol());
      Rcpp::List result = solve_QP(Rcpp::Named("Dmat") = D_nng_R,
                                   Rcpp::Named("dvec") = d_R2,
                                   Rcpp::Named("Amat") = Amat_R,
                                   Rcpp::Named("bvec") = b0);
      Rcpp::NumericVector coef_nng_R = result["solution"];
      best_coef_nng  = Rcpp::as<arma::vec>(coef_nng_R);
    }

    beta_nng_arma = beta_arma % best_coef_nng;
    Rcpp::NumericVector beta_nng = Rcpp::wrap(beta_nng_arma);

    return beta_nng;
  }

  
private:
  Rcpp::List h_j_list;
  int p;
  Rcpp::List rho_list;
  Rcpp::IntegerVector mi;
  std::vector<arma::mat> Psi_j_list;
  double tau0_sigma0;
};


//////////////////////////////// Export BETA
// Define a global pointer to an BETA instance
Rcpp::XPtr<BETA>* BETA_instance_ptr = nullptr;

// Function to initialize the BETA instance
// [[Rcpp::export]]
void initialize_BETA_instance(Rcpp::List h_j_list, int p, Rcpp::List rho_list, Rcpp::IntegerVector mi) {
  BETA* ptr = new BETA(h_j_list, p, rho_list, mi);
  BETA_instance_ptr = new Rcpp::XPtr<BETA>(ptr, true);
} 

// Export r_j_cpp
// [[Rcpp::export]]
Rcpp::List r_j_cpp_R(Rcpp::List U_j_list, Rcpp::IntegerVector me_num) {
  if (BETA_instance_ptr == nullptr) {
    Rcpp::stop("BETA instance is not initialized. Call initialize_BETA_instance first.");
  }
  return BETA_instance_ptr->get()->r_j_cpp(U_j_list, me_num);
}

// Export beta_ele_cpp
// [[Rcpp::export]]
Rcpp::List beta_ele_cpp_R(const Rcpp::NumericMatrix& U, const Rcpp::NumericVector& R,
                          double lambda, int replicate, int n,
                          const Rcpp::NumericVector& y) {
  if (BETA_instance_ptr == nullptr) {
    Rcpp::stop("BETA instance is not initialized. Call initialize_BETA_instance first.");
  }
  return BETA_instance_ptr->get()->beta_ele_cpp(U, R, lambda, replicate, n, y);
}

// Export beta_cpp
// [[Rcpp::export]]
Rcpp::NumericVector beta_nng_cpp_R(const Rcpp::List& beta_ele, int replicate, int n,
                                   const Rcpp::NumericVector& y, const Rcpp::NumericMatrix& Amat,
                                   double s2) {
  if (BETA_instance_ptr == nullptr) {
    Rcpp::stop("BETA instance is not initialized. Call initialize_BETA_instance first.");
  }
  return BETA_instance_ptr->get()->beta_nng_cpp(beta_ele, replicate, n, y, Amat, s2);
}



